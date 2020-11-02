from itertools import permutations
from typing import AnyStr, Optional
from fastapi import FastAPI, File, UploadFile, HTTPException
from .utils.file_to_str import file_to_str
from .utils.is_valid import is_valid_dna
from .rosalind.recv import reverse_complement as rc
from .utils.read_fasta import fasta_to_sequence
from .utils.str_to_nums import str_to_nums

from .routers import problems

from itertools import product

app = FastAPI(
    title="Rosalind Answer API",
    description="An API that implements my solutions to the Rosalind problem set"
)

app.include_router(problems.router, prefix="/problems")

from math import factorial


def read_root():
    return {"error": "Endpoint not implemented yet"}


