from fastapi import FastAPI
from starlette.requests import Request
from starlette.responses import RedirectResponse

from .routers import problems

app = FastAPI(
    title="Franklin API",
    description="""My solutions to the [Rosalind](http://rosalind.info/) problem set, grouped by tag (some endpoints have multiple tags). 
    View the [source code here](https://github.com/boonepeter/franklin-api) (Currently private to not distribute solutions)."""
)

app.include_router(problems.router, prefix="/problems")

@app.get("/")
async def redirect():
    """
    Redirect requests for '/' to the docs
    """
    response = RedirectResponse(url="/docs")
    return response



@app.middleware("http")
async def case_sens_middleware(request: Request, call_next):
    path = request.scope["path"].lower()
    request.scope["path"] = path
    response = await call_next(request)
    return response
