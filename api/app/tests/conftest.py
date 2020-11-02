import os

def pytest_addoption(parser):
    parser.addoption("--route", action="store", default="default_name")

def pytest_generate_tests(metafunc):
    dirpath = os.path.join(os.path.dirname(__file__), "data")
    files = [f for f in os.listdir(dirpath)]
    routes = set()
    for f in files:
        name = f.split("_")[1]
        routes.add(name)

    option_value = metafunc.config.option.route
    if "route" in metafunc.fixturenames and routes is not None:
        metafunc.parametrize("route", routes)
