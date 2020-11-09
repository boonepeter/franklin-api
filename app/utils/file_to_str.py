from typing import OrderedDict, List
from fastapi import File, UploadFile
from fastapi.exceptions import HTTPException
from .read_fasta import fasta_to_sequence
import os

async def file_to_str(file: UploadFile = File(...), max_size: int=10000) -> str:
    if file.__sizeof__() > max_size:
        print(file.__sizeof__())
    contents = await file.read()
    contents = contents.decode("utf-8").strip()
    #if len(contents) == 0:
    #    raise HTTPException(status_code=500, detail="File is empty")
    return contents


async def file_to_fastas(file: UploadFile=File(...), expected: int=-1, max_size: int=10000) -> OrderedDict[str, str]:
    contents = await file_to_str(file, max_size)
    seqs = fasta_to_sequence(contents)
    if not expected == -1 and not len(seqs) == expected:
        raise HTTPException(status_code=500, detail=f"Expected {expected} sequence(s), received {len(seqs)}")
    return seqs


async def file_to_seqs(file: UploadFile=File(...), expected: int=-1, max_size: int=10000) -> List[str]:
    fastas = await file_to_fastas(file, expected=expected, max_size=max_size)
    return list(fastas.values())