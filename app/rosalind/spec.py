

from typing import Dict, List
from .prtm import masses

def calc_protein(prefix_masses: List[float]) -> str:
    mass_to_protein: Dict[float, str] = {}
    for k, v in masses.items():
        mass_to_protein[v] = k
    build = ""
    for i in range(1, len(prefix_masses)):
        w = prefix_masses[i] - prefix_masses[i - 1]
        key = min(mass_to_protein, key=lambda x: abs(x - w))
        build += mass_to_protein[key]
    return build

