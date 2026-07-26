"""
Filename: mosfet_selection.py

Description:
    Calculates a dimensional quality score for 
    selecting MOSFETs based on their parameters 
    for the triple bridge driver.

    NOTE:
    MOSFET costs are in AUD, but the script will work 
    with any currency as long as it is used consistently.
    
    Also they must all be either n-channel or p-channel or
    LGBT. Cannot mix types.
    
    NOTE:
    Selected:
    Mosfet("SUD90330E", 3.4 * nullset, 32 * nC, 37.5 * mR, 35.1 * current, 200 * voltage)
"""

from __future__ import annotations

from dataclasses import dataclass
from operator import itemgetter

from picounits import Q
from picounits.constants import nullset, charge, resistance, current, voltage, nano, milli


# Defines the mosfet set for usage
component_set: list[Mosfet] = []

@dataclass(slots=True, frozen=True)
class Mosfet:
    """ Stores the values for a specific mosfet """
    name: str
    cost: Q
    gate_charge: Q          # Maximum total gate charge
    on_resistance: Q        # Static drain to source on resistance
    current: Q              # Maximum drain current
    max_voltage: Q          # Maximum source to drain voltage

    def __post_init__(self) -> None:
        """ Appends self to the mosfet set """
        component_set.append(self)

    def score(self) -> Q:
        """ Calculates a floating point score for the component """
        denom = (6 * self.cost) * self.gate_charge * self.on_resistance
        value = self.current * self.max_voltage / denom

        # Returns the resulting value with units
        return  value

    @property
    def _name(self):
        """ Returns name as attributes """
        return (
            f"<{self.name}($={6 * self.cost:.2f}, Q_g={self.gate_charge:.2f}, "
            f"R_on={self.on_resistance:.2f}, "
            f"I_373k={self.current:.2f}, v_m={self.max_voltage:.2f})>"
        )

    def __repr__(self) -> str: return self._name
    def __str__(self) -> str: return self._name


""" Define your MOSFETs here. """
nC = 1 * nano * charge
mR = 1 * milli * resistance

# All are n-channel mosfet using DPAK Style package (most values are at v_gs @ 10V, 25C, cut-tape)

# These below are from mouser catalog
Mosfet("FDD18N20LZ", 4.41 * nullset, 40 * nC, 125 * mR, 16 * current, 200 * voltage)
Mosfet("FQD18N20V2", 3.18 * nullset, 26 * nC, 140 * mR, 15 * current, 200 * voltage)
Mosfet("RCJ200N20",  3.57 * nullset, 40 * nC, 130 * mR, 20 * current, 200 * voltage)
Mosfet("RD3T100CN", 3.58 * nullset, 25 * nC, 182 * mR, 10 * current, 200 * voltage)
Mosfet("RD3R02BBH", 3.74 * nullset, 12.4 * nC, 81 * mR, 20 * current, 150 * voltage)
Mosfet("fdd86250", 3.22 * nullset, 33 * nC, 22 * mR, 8 * current, 150 * voltage)
Mosfet("IRFR24N15DPbF", 2.55 * nullset, 45 * nC, 95 * mR, 24 * current, 150 * voltage)
Mosfet("SQD25N15-52", 7.45 * nullset, 51 * nC, 52 * mR, 25 * current, 150 * voltage)
Mosfet("DMN15H310SK3", 1.52 * nullset, 4.6 * nC, 310 * mR, 8.3 * current, 150 * voltage)
Mosfet("IPD110N12N3", 3.92 * nullset, 65 * nC, 11 * mR, 75 * current, 120 * voltage)

# These below are from digikey catalog
Mosfet("IRF640NPbF", 3.17 * nullset, 67 * nC, 150 * mR, 18 * current, 200 * voltage)
Mosfet("IRFS4620TRLPBF", 5.2 * nullset, 38 * nC, 77.5 * mR, 24 * current, 200 * voltage)
Mosfet("SUM90142E-GE3", 6.8 * nullset, 87 * nC, 16 * mR, 90 * current, 200 * voltage)
Mosfet("FQD7N20L", 2.10 * nullset, 9 * nC, 750 * mR, 5.5 * current, 200 * voltage)
Mosfet("SUD90330E", 3.4 * nullset, 32 * nC, 37.5 * mR, 35.1 * current, 200 * voltage)
Mosfet("IRFS4227TRLPBF", 6.84 * nullset, 98 * nC, 26 * mR, 62 * current, 200 * voltage)

# These below are from element14 catalog
Mosfet("IRFR220NPbF", 1.370 * nullset, 23 * nC, 600 * mR, 5 * current, 200 * voltage)


# Voltage too low for 96 V operation (Mouser)
# Mosfet("IRFS4010TRLPBF", 8.13 * nullset, 215 * nC, 4.7 * mR, 180 * current, 100 * voltage)

""" Displays the results as a "score" table """

results: list[tuple[Mosfet, Q]] = []

for component in component_set:
    score = component.score()
    results.append((score, component))

# Sorts by score (highest to lowest) & prints
results.sort(key=itemgetter(0), reverse=True)

print("--" * 30, " Score Table & Components (H2L) ", "--" * 30)
for result in results:
    # Displays the, score and component
    print(f"Score: {result[0]}, Component: {result[1]}")
