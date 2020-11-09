import os

from ..main import app

def pytest_generate_tests(metafunc):
    dirpath = os.path.join(os.path.dirname(__file__), "data")
    files = [f for f in os.listdir(dirpath)]
    file_prefixes = set()
    for f in files:
        split = f.split("_")
        name = split[1] + "_" + split[2]
        file_prefixes.add(name)

    skip = ["/openapi.json", "/docs", "/docs/oauth2-redirect", "/redoc", "/"]

    routes = [i.path for i in app.routes if i.path not in skip]

    if "path" in metafunc.fixturenames and routes is not None:
        metafunc.parametrize("path", routes)

    if "file_prefix" in metafunc.fixturenames and file_prefixes is not None:
        metafunc.parametrize("file_prefix", file_prefixes)
