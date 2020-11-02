from fastapi import FastAPI, File, UploadFile
import os
async def file_to_str(file: UploadFile = File(...), max_size: int=10000) -> str:
    if file.__sizeof__() > max_size:
        print(file.__sizeof__())
    contents = await file.read()
    return contents.decode("utf-8").strip()
