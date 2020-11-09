
from typing import List
from fastapi import HTTPException

def assert_len(expected: int, actual: List, error: bool=True) -> bool:
    if not expected == len(actual):
        if error:
            raise HTTPException(status_code=500, detail=f"Expected {expected} argument(s), received {len(actual)}")
        return False
    else:
        return True