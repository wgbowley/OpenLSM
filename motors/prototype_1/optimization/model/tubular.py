"""
Filename: tubular.py

Description:
    Thermal and magnetic model of a tubular linear motor
    for usage in co-simulation. 
    
    NOTE:
    Slots parameters define pole spacing during optimization.
    
    NOTE:
    Thermal model isn't intended to be usage in mechanical-thermal
    model instead static position transient thermal should be used. 
"""

from __future__ import annotations

from math import ceil
from pathlib import Path
from pyfea import ampere, meter, millimeter, watt, dimensionless, Quantity as Q

from pyfea.domain.units import Parser, Configuration
from pyfea.domain.materials.manager import MaterialManager

from pyfea.domain.geometry.definitions import CoordinateSystem
from pyfea.domain.geometry.builder import Builder, VectorGeometry
from pyfea.domain.geometry.elements.vectors import CSGNode
from pyfea.domain.geometry.domain import Domain, BoundaryType

from pyfea.solver.solver_interface import BaseSolver, MagneticSolver, ThermalSolver
from pyfea.domain.geometry.elements.metadata import MagneticData, ThermalData
from pyfea.domain.circuits.builder import Circuit, Configuration as CircuitConfig


class ModelError(Exception):
    """ Exception for tubular motor modelling errors """
    def __init__(self, func_name: str = "", error: str = "") -> None:
        """ Returns a error message with caller and the error """
        msg = f"{func_name} raised error: {error}"
        super().__init__(msg)


class TubularLinearMotor:
    """ CSG model for a tubular motor with domain construction method """
    # Simulation element tags
    ENVIRONMENT_ID = 0 * dimensionless
    SLOT_ID =   1 * dimensionless
    CORE_ID = 2 * dimensionless
    POLE_ID = 3 * dimensionless
    TUBE_ID = 4 * dimensionless 
    
    # Phases within the motor
    PHASES = [
        Circuit("phase_a", 0 * ampere, CircuitConfig.SERIES), 
        Circuit("phase_b", 0 * ampere, CircuitConfig.SERIES),
        Circuit("phase_c", 0 * ampere, CircuitConfig.SERIES)
    ]
    
    def __init__(self, configuration_path: Path) -> None:
        """ Initializes the class & defines dependencies """
        self.params = self._get_parameters(configuration_path)

        # Defines the coordinate system and model materials
        self.coordinate_system = CoordinateSystem.AXI_SYMMETRIC
        self._load_material()
        
        # Initializes variables for future usage
        self.total_poles = 0
        self.slot_pitch = 0.0
        self.pole_pitch = 0.0
        self.armature_length = 0.0
        self.effective_length = 0.0
    
    def construct_domain(self, solver: BaseSolver) -> Domain:
        """ Constructs the domain based on solver physics domain """
        solver_interfaces = solver.__class__.__bases__
    
        # Derived parameters & reloads materials
        self._derived_parameters()
        self._load_material()
    
        for solver in solver_interfaces:
            if solver == ThermalSolver:
                return ConstructThermal.build(self)
            
            if solver == MagneticSolver:
                return ConstructMagnetic.build(self)
            
        msg = f"{solver_interfaces!r} is not supported by {self.__class__.__name__}"
        raise ModelError("TubularLinearMotor.construct_domain()", msg)
        
    def update_parameters(
        self, 
        pole_slot_ratio: Q, 
        pole_grade: str, 
        wire_diameter: Q, 
        slot_axial_length: Q,
        slot_outer_radius: Q, 
        slot_axial_spacing: Q
    ) -> None:
        """ Updates parameters within the configuration to reflect changes """
        self.params.find_and_replace("model.number_pairs", pole_slot_ratio[0])
        self.params.find_and_replace("model.number_slots", pole_slot_ratio[1])
        self.params.find_and_replace("stator_poles.grade", pole_grade)
        self.params.find_and_replace("armature_slots.wire_diameter", wire_diameter)
        self.params.find_and_replace("armature_slots.axial_length", slot_axial_length)
        self.params.find_and_replace("armature_slots.outer_radius", slot_outer_radius)
        self.params.find_and_replace("armature_core.outer_radius", slot_outer_radius)
        self.params.find_and_replace("armature_core.axial_slot_spacing", slot_axial_spacing)
        
    def build_armature(self) -> tuple[CSGNode, list[VectorGeometry]]:
        """ Builds the armature via constructing the core and its slots """
        parameters = self.params
        core_radial_thickness = (
            parameters.armature_core.outer_radius - parameters.armature_core.inner_radius
        )
        
        # Constructs the rectangular base shape to subtract the slots against.
        core = Builder.create_rectangle(
            (
                parameters.armature_core.inner_radius,
                - self.armature_length / 2
            ),
            core_radial_thickness, self.armature_length
        )
        
        # Constructs the slots and remove the slot area from the armature
        slots = []
        for slot in range(0, parameters.model.number_slots.value):
            offset = - self.effective_length / 2
            bottom_left = offset + slot * self.slot_pitch
            
            slot = Builder.create_rectangle(
                (parameters.armature_slots.inner_radius, bottom_left),
                parameters.armature_slots.outer_radius - parameters.armature_slots.inner_radius,
                parameters.armature_slots.axial_length
            )
            core.subtract(slot)
            slots.append(slot)
            
        return core, slots
    
    def build_stator(self, number_poles: Q) -> tuple[list, VectorGeometry]:
        """ Builds the stator poles and their enclosing tube """
        poles = []
        for pole in range(0, number_poles.value):
            offset = - number_poles * self.pole_pitch / 2
            bottom_left = offset + pole * self.pole_pitch
            
            pole_shape = Builder.create_rectangle(
                (0 * millimeter, bottom_left),
                self.params.stator_poles.outer_radius, self.pole_pitch
            )
            poles.append(pole_shape)

        tube = Builder.create_rectangle(
            (self.params.stator_poles.outer_radius, - number_poles * self.pole_pitch / 2),
            self.params.stator_tube.outer_radius - self.params.stator_tube.inner_radius,
            number_poles * self.pole_pitch
        )
        return poles, tube

    def build_boundary(self, number_poles: Q) -> VectorGeometry:
        """ Builds the boundary shape via parameters with 20% margin for safety """
        return Builder.create_rectangle(
            (0 * millimeter, - 1.2 * number_poles * self.pole_pitch / 2),
            self.params.armature_core.outer_radius * 1.2,
            number_poles * self.pole_pitch * 1.2
        )
    
    def _get_parameters(self, path: Path) -> Parser:
        """ Attempts to resolve user configuration file. Default to model path """
        target = Path(path).expanduser().resolve()
        
        if not target.exists():
            msg = f"{path} does not exist or could not be resolved."
            raise ModelError("TubularLinearMotor._get_parameters()", msg)

        return Parser.open(target, Configuration)
    
    def _load_material(self) -> None:
        """ Builds out the material manager with materials """
        manager = MaterialManager()
        
        # Finds the material in the .uiv material library
        self.environmental_material = manager.use_material(
            self.params.model.environmental_material
        )
        self.armature_core_material = manager.use_material(
            self.params.armature_core.material
        )
        self.armature_slots_material = manager.use_material(
            self.params.armature_slots.material
        )
        self.stator_tube_material = manager.use_material(
            self.params.stator_tube.material
        )
        self.stator_poles_material = manager.use_material(
            self.params.stator_poles.material, 
            grade = self.params.stator_poles.grade
        )
        
    def _derived_parameters(self) -> None:
        """ Calculates derived parameters from base parameters """
        # Segment (number under the armature), boundary (approx for infinite track)
        segment_poles = 2 * self.params.model.number_pairs
        boundary_poles = 4 * self.params.model.boundary_pairs

        self.total_poles = segment_poles + boundary_poles
        
        # Slot pitch is the distance between adjacent slots (start 1 -> start 2)
        axial_spacing = self.params.armature_core.axial_slot_spacing
        self.slot_pitch = self.params.armature_slots.axial_length + axial_spacing
        
        # Physical length and effective magnetic length of the armature
        self.armature_length = self.slot_pitch * self.params.model.number_slots
        self.effective_length = self.armature_length - axial_spacing
        
        # Adds another axial_spacing to the armature length for equivalents end-segments
        self.armature_length += axial_spacing
        
        # Pole pitch is the distance between adjacent poles (start 1 -> start 2)
        self.pole_pitch = self.effective_length / segment_poles
        
        pole_length = self.params.stator_poles.axial_length
        if self.pole_pitch < pole_length:
            msg = "Failed to derive parameters, overlapping stator poles: "
            msg += f"{self.pole_pitch:.3f} : {pole_length:.3f}"
            raise ModelError("TubularLinearMotor._derived_parameters()", msg)
        

class ConstructMagnetic:
    """ Constructs magnetic domain via adding magnetic metadata to geometry """
    @classmethod
    def calculate_number_turns(cls, motor: TubularLinearMotor) -> int:
        """ Calculates the approximate number of turns within a slot """
        parameters = motor.params
        radius = parameters.armature_slots.outer_radius - parameters.armature_slots.inner_radius

        # Calculates the cross sectional area, wire area and than effective material area
        slot_area = radius * parameters.armature_slots.axial_length
        wire_area = parameters.armature_slots.wire_diameter ** 2

        effective_area = slot_area * parameters.armature_slots.fill_factor
        turns = ceil(effective_area / wire_area)
        
        if turns < 0:
            msg = f"Derived parameter 'turns' cannot be {turns}. Slots must have non-zero area"
            raise ModelError("ConstructMagnetic.calculate_number_turns", msg)
        
        return turns
        
    @classmethod
    def build(cls, motor: TubularLinearMotor) -> Domain:
        """ Builds the magnetic simulation via adding metadata to geometry """
        parameters = motor.params
        parts = []
        
        # Builds armature and stator than extract parts
        core, slots = motor.build_armature()
        poles, tube = motor.build_stator(motor.total_poles)
        boundary    = motor.build_boundary(motor.total_poles)
        
        # Defines simulation parts via promoting and metadata
        parts = []

        # Adds the slots to the domain
        turns = cls.calculate_number_turns(motor)
        for index, slot in enumerate(slots):
            # Sets phase of slot in pattern [a, b, c] with alternating polarity for slots
            phase = motor.PHASES[index % len(motor.PHASES)]
            polarity = +1 if index % 2 == 0 else -1
            
            # Constructs meta-data and promotes to part while appending to domain
            meta = MagneticData(
                motor.SLOT_ID, motor.armature_slots_material, phase, turns * polarity,
                parameters.armature_slots.wire_diameter
            )
            parts.append(Builder.promote_to_part(slot, meta))

        # Adds the poles to the domain
        for index, pole in enumerate(poles):
            # Alternate magnetization direction every pole (e.g., N-S-N-S)
            pole_magnetization = 90 if index % 2 == 0 else - 90

            # Constructs meta-data and promotes to part while appending to domain
            meta = MagneticData(
                motor.POLE_ID, motor.stator_poles_material,
                magnetization = pole_magnetization * dimensionless
            )
            parts.append(Builder.promote_to_part(pole, meta))

        # Adds the armature core and the stator tube to the domain
        parts.append(
            Builder.promote_to_part(
                core, MagneticData(motor.CORE_ID, motor.armature_core_material)
            )
        )
        parts.append(
            Builder.promote_to_part(
                tube, MagneticData(motor.TUBE_ID, motor.stator_tube_material)
            )
        )
        
        # Overall simulation problem definition
        meta = MagneticData(motor.ENVIRONMENT_ID, motor.environmental_material)
        return Domain(
            parts, BoundaryType.DIRICHLET, meta, motor.coordinate_system, boundary
        )


class ConstructThermal:
    """ Constructs thermal domain via adding thermal metadata to geometry """
    @classmethod
    def build(cls, motor: TubularLinearMotor) -> Domain:
        """ Builds the thermal simulation via adding metadata to geometry """
        parameters = motor.params
        parts = []
        
        # Builds armature and stator than extract parts
        core, slots = motor.build_armature()
        poles, tube = motor.build_stator(motor.params.model.number_pairs * 2)
        boundary    = motor.build_boundary(motor.params.model.number_pairs * 2)
        
        # Defines simulation parts via promoting and metadata
        parts = []
        
        # Adds the slots to the domain
        meta = ThermalData(
            motor.SLOT_ID, motor.armature_slots_material,
            volumetric_heating = 0 * watt / meter ** 3
        )
        for slot in slots: 
            parts.append(Builder.promote_to_part(slot, meta))
        
        # Adds the poles to the domain
        meta = ThermalData(motor.POLE_ID, motor.stator_poles_material)
        for pole in poles: 
            parts.append(Builder.promote_to_part(pole, meta))
            
        # Adds the armature core and the stator tube to the domain
        parts.append(
            Builder.promote_to_part(
                core, ThermalData(motor.CORE_ID, motor.armature_core_material)
            )
        )

        parts.append(
            Builder.promote_to_part(
                tube, ThermalData(motor.TUBE_ID, motor.stator_tube_material)
            )
        )
        
        # Overall simulation problem definition
        meta = ThermalData(
            motor.ENVIRONMENT_ID, motor.environmental_material,
            temperature=motor.params.thermal.atmospheric_temperature,
            convection_coefficient=motor.params.thermal.convection_coefficient
        )
        return Domain(
            parts, BoundaryType.CONVECTION, meta, motor.coordinate_system, boundary
        )
        