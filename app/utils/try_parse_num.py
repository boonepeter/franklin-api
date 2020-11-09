from fastapi import HTTPException

def try_parse_float(in_string: str) -> float:
    try:
        x = float(in_string)
        return x
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to parse float: {in_string}, {e}")

def try_parse_int(in_string: str) -> int:
    try:
        x = int(in_string)
        return x
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to parse int: {in_string}, {e}")