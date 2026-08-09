"""
Filename: model.py

Description:
    Magnet model of an tubular ironless pcb armature linear
    motor for usage in FEM simulation
"""

from pyfea import nullset, ampere
from pyfea.domain.units import DynamicLoader

from pyfea.domain.materials.manager import MaterialManager
from pyfea.domain.geometry.definitions import CoordinateSystem

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
        self.pole_pitch = self.params.stator.dipole.axi_length
        self.number_poles = self.params.stator.tube.axi_length / self.pole_pitch

        # Computes the slot pitch from pole pitch
        self.slot_pitch = self.pole_pitch / 3

        # Computes the number of pcb's within each slot
        self.pcb_layers = self.params.armature.sub_slot.axi_layers
        self.pcb_trace_thickness = self.params.armature.sub_slot.trace_thickness
        self.pcb_substrate = self.params.armature.sub_slot.axi_len_substrate

        self.pcb_axi_thickness = self.pcb_layers * self.pcb_trace_thickness + self.pcb_substrate
        self.sub_slots = self.slot_pitch // self.pcb_axi_thickness

        # Computes the armature length
        self.armature_length = self.slot_pitch * self.params.armature.core.number_slots

        # Computes the radial placement of components
        self.pole_radius = self.params.stator.dipole.rad_thickness
        self.tube_radius = self.pole_radius + self.params.stator.tube.rad_thickness
        self.arm_inner_radius = self.tube_radius + self.params.armature.core.rad_clearance
        self.winding_inner_radius = self.arm_inner_radius + self.params.armature.core.pcb_edge
        self.winding_outer_radius = self.winding_inner_radius + self.params.armature.slot.rad_thickness
