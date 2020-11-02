from typing import List
from fastapi import HTTPException


def str_to_nums(input_str: str, split_char: str=None, expected: int=-1) -> List[int]:
    split = input_str.split(split_char)
    if not expected == -1 and not len(split) == expected:
        raise HTTPException(status_code=500, detail=f"Expected {expected} arguments, received {len(split)}")
    try:
        nums = [int(i) for i in split]
    except ValueError as e:
        raise HTTPException(status_code=500, detail=f"Invalid argument, expected integers: {e}")
    return nums
