import re
from typing import Match

def float_replace(match: Match):
    s = match.group(0)
    f = round(float(s), 4)
    replaced = f"{f:.3f}"
    return replaced

def replace_floats(in_str: str) -> str:
    f = r"[+-]?(?=\d*[.eE])(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?"
    replaced = re.sub(f, float_replace, in_str)
    return replaced