
from fastapi import HTTPException
from typing import List


def split_str(input_string: str, split_char: str=None, expected: int=-1) -> List[str]:
    split_list = input_string.split(split_char)
    if not expected == -1:
        if not len(split_list) == expected:
            raise HTTPException(status_code=500, detail=f"Expected {expected} argument(s), received {len(split_list)}")
    return split_list
