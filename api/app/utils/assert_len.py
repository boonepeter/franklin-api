
from fastapi import HTTPException

def assert_len(expected: int, actual: int, error: bool=True) -> bool:
    if not expected == actual:
        if error:
            raise HTTPException(status_code=500, detail=f"Expected {expected} argument(s), received {actual}")
        return False
    else:
        return True