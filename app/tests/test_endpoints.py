from fastapi import UploadFile, File
from fastapi.testclient import TestClient
import os
import io

from starlette.routing import BaseRoute
from ..main import app

from ..routers.problems import Answer, router
from .util import replace_floats

"""
def test_no_file(path: str):
    with TestClient(app) as client:
        response = client.post(path)
        assert response.status_code == 422    
    return

def test_empty_file(path: str):
    with TestClient(app) as client:
        file = io.StringIO("xyz")
        files = {"file": file}
        response = client.post(path, files=files)
        assert response.status_code == 500
"""
def test_docs(path: str):
    r = [r for r in app.routes if r.path == path]
    route = r[0]
    assert route is not None
    assert hasattr(route, "description")
    assert hasattr(route, "path")
    if hasattr(route, "description") and hasattr(route, "path"):
        description: str = route.description
        splits = description.split("```")
        
        if len(splits) <= 2:
            print(f"{path} not implemented")
            return
        if len(splits) > 2:
            dataset = splits[1].strip()
            output = splits[3].strip()
            file = io.StringIO(dataset)
            files = { 'file': file }
            with TestClient(app) as client:
                response = client.post(path, files=files)
                assert response.status_code == 200
                json = response.json()
                if "answer" not in json:
                    assert json["error"] == "not implemented"
                    return
                assert "answer" in json
                if json["answer"] == "":
                    print(f"{path} not implemented")
                    return
                replaced = replace_floats(str(json['answer']))
                r_output = replace_floats(output)
                assert replaced == r_output



def test_datasets(file_prefix: str):
    with TestClient(app) as client:
        data_folder = os.path.join(os.path.dirname(__file__), "data")
        f = os.path.join(data_folder, f"rosalind_{file_prefix}_dataset.txt")
        out = os.path.join(data_folder, f"rosalind_{file_prefix}_output.txt")
        name = file_prefix.split("_")[0]
        url = f"/problems/{name.upper()}"
        with open(f, 'rb') as tmp:
            files = {'file': tmp}
            response = client.post(url, files=files)
            assert response.status_code == 200
            json = response.json()
            if "answer" not in json:
                assert json["error"] == "not implemented"
                return
            assert "answer" in json
            with open(out, "r") as out_file:
                replaced = replace_floats(str(json["answer"]))
                replaced_output = replace_floats(out_file.read())
                assert replaced == replaced_output
