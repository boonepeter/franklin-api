from fastapi import UploadFile, File
from fastapi.testclient import TestClient
import os

from ..main import app

client = TestClient(app)


def test_problem(route: str):
    with TestClient(app) as client:
        data_folder = os.path.join(os.path.dirname(__file__), "data")
        f = os.path.join(data_folder, f"rosalind_{route}_1_dataset.txt")
        out = os.path.join(data_folder, f"rosalind_{route}_1_output.txt")
        url = f"/problems/{route.upper()}"
        with open(f, 'rb') as tmp:
            files = {'file': tmp}
            response = client.post(url, files=files)
            assert response.status_code == 200
            json = response.json()
            if "answer" in json:
                assert "answer" in json
                with open(out, "r") as out_file:
                    assert str(json["answer"]) == out_file.read()
            else:
                assert "error" in json


