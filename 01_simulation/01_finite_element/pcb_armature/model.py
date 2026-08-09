"""
Filename: model.py

Description:
    Magnet model of an tubular ironless pcb armature linear
    motor for usage in FEM simulation
"""

from pyfea import nullset, ampere, mm
from pyfea.domain.units import DynamicLoader

from pyfea.domain.materials.manager import MaterialManager
from pyfea.domain.geometry.definitions import CoordinateSystem
from pyfea.domain.geometry.builder import Builder, VectorGeometry, Part
from pyfea.domain.geometry.domain import Domain, BoundaryType

from pyfea.solver.solver_interface import BaseSolver, MagneticSolver
from pyfea.domain.circuits.builder import StaticCircuit, Configuration as CircuitConfig


class ModelError(Exception):
    """ Exception for tubular motor modelling errors """
    def __init__(self, func_name: str = "", error: str = "") -> None:
        """ Returns a error message with caller and the error """
        msg = f"{func_name} raised error: {error}"
        super().__init__(msg)


class PCBArmatureMotor:
    """ CSG Model for a pcb armature linear motor """
    ENVIRONMENT_ID = 0 * nullset
    STATOR_ID = 1 * nullset
    ARMATURE_ID = 2 * nullset

    # Static independent circuits non-connected to a circuit solver
    PA = StaticCircuit("PA", 0 * ampere, CircuitConfig.series)
    PB = StaticCircuit("PB", 0 * ampere, CircuitConfig.series)
    PC = StaticCircuit("PC", 0 * ampere, CircuitConfig.series)

    def __init__(self, parameters: DynamicLoader) -> None:
        """ Initializes the class & defines dependencies  """
        self.params = parameters
        self.coordinate_system = CoordinateSystem.AXI_SYMMETRIC

        # Loads the material properties and computes derived parameters
        self._load_materials_properties()
        self._compute_derived_parameters()

    def construct_domain(self, solver: BaseSolver) -> tuple[Domain, list[Part]]:
        """ Constructs the domain based on solver physics """
        solver_interfaces = solver.__class__.__bases__

        for solver in solver_interfaces:
            if solver == MagneticSolver:
                return None # Temp.

        msg = f"{solver_interfaces!r} is not supported by {self.__class__.__name__!r}"
        raise ModelError("PCBArmatureMotor.construct_domain()", msg)

    def build_armature(self) -> list[VectorGeometry]:
        """ Builds the slots within the armature """
        slots = []

        number_layers = self.sub_slots * self.number_slots * self.pcb_layers
        for index in range(0, int(number_layers.value)):
            offset = - self.armature_length / 2
            bottom_left = offset + index * (self.pcb_trace_thickness + self.pcb_substrate)

            slot_shape = Builder.rectangle(
                (self.winding_inner_radius, bottom_left),
                self.winding_thickness, self.pcb_trace_thickness
            )
            slots.append(slot_shape)

        return slots

    def build_stator(self) -> list[VectorGeometry]:
        """ Builds the poles within the stator """
        dipoles = []

        for index in range(0, int(self.number_dipoles.value)):
            offset = - self.stator_axi_len / 2
            bottom_left = offset + index * self.dipole_pitch + self.dipole_pitch / 2

            dipoles_shape = Builder.rectangle(
                (0 * mm, bottom_left), self.dipole_radius, self.dipole_pitch
            )
            dipoles.append(dipoles_shape)

        return dipoles

    def build_boundary(self) -> VectorGeometry:
        """ Builds the boundary shape via tube length with 20% margin axially and 100% radially """
        max_length = 1.2 * self.params.stator.tube.axi_length
        max_radius = 2 * self.winding_outer_radius

        return Builder.rectangle((0 * mm, -max_length / 2), max_radius, max_length)

    def _load_materials_properties(self) -> None:
        """ Loads the material properties into class data """
        manager = MaterialManager()

        # Finds the materials in the .uiv material library
        self.env_material = manager.use_material(self.params.model.environmental_material)
        self.stator_material = manager.use_material(self.params.armature.sub_slot.material)

        self.armature_material = manager.use_material(
            self.params.stator.dipole.material,
            grade=self.params.stator.dipole.grade
        )

    def _compute_derived_parameters(self) -> None:
        """ Computes derived parameters from static parameters """
        self.stator_axi_len = self.params.stator.tube.axi_length
        self.dipole_pitch = self.params.stator.dipole.axi_length

        # Computes the number of dipoles within the stator
        self.number_dipoles = self.stator_axi_len / self.dipole_pitch

        # Computes the slot pitch from pole pitch
        self.slot_pitch = self.dipole_pitch / 3

        # Computes the number of pcb's within each slot
        self.pcb_layers = self.params.armature.sub_slot.axi_layers
        self.pcb_trace_thickness = self.params.armature.sub_slot.trace_thickness
        self.pcb_substrate = self.params.armature.sub_slot.axi_len_substrate

        self.pcb_axi_thickness = self.pcb_layers * self.pcb_trace_thickness + self.pcb_substrate
        self.sub_slots = self.slot_pitch // self.pcb_axi_thickness

        # Computes the armature length
        self.number_slots = self.params.armature.core.number_slots
        self.armature_length = self.slot_pitch * self.number_slots

        # Computes the radial placement of components
        self.dipole_radius = self.params.stator.dipole.rad_thickness
        self.tube_radius = self.dipole_radius + self.params.stator.tube.rad_thickness
        self.arm_inner_radius = self.tube_radius + self.params.armature.core.rad_clearance
        self.winding_inner_radius = self.arm_inner_radius + self.params.armature.core.pcb_edge

        self.winding_thickness = self.params.armature.slot.rad_thickness
        self.winding_outer_radius = self.winding_inner_radius + self.winding_thickness
