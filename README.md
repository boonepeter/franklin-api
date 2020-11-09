# Franklin API

This [FastAPI](https://fastapi.tiangolo.com/) API provides the answers to problems in the [Rosalind](http://rosalind.info/) problem set.

A deployed version of this API can be found [here hosted on Azure](https://boone-rosalind.azurewebsites.net/docs). You can view the docs and test out some of the endpoints there.

__Don't cheat!__

> I set this up as a personal project to explore FastAPI and solve some of the Rosalind problems. While you could use this to cheat on the Rosalind problems, I urge you not to. The whole point of the Rosalind problem set is to learn the algorithms and concepts behind some bioinformatics problems.

## Endpoints

All endpoints accept `POST` requests with `multipart/form-data` as the body. This data should be just one [`text/plain`](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types) file. You can hit this endpoint using `requests` like this:

```python
import requests

local_filepath = "./path/to/local/file.txt"
url = "https://boone-rosalind.azurewebsites.net/problems/dna"

with open(local_filepath, 'rb') as file:
    files = {"file": file}
    response = requests.post(url, files=files)
    print(response.json())
```

If you have the dataset as a string, you can use `StringIO` to create a file-like object:

```python
import requests
import io

contents = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
url = "https://boone-rosalind.azurewebsites.net/problems/dna"
file = io.StringIO(contents)
files = {"file": file}
response = requests.post(url, files=files)
print(response.json())
```

You can see the format of the file that each endpoint accepts in each endpoint description. The formats match the datasets that Rosalind provides. With a `200` status code, each endpoint will respond with the following json:

```json
// a string
{
    "answer": "ACGT"
}

// or a number

{
    "answer": 0.5
}
```

You can write the contents of `response.json()["answer"]` to a file, and that should be the answer that Rosalind expects.

## Testing

![Test](https://github.com/boonepeter/franklin-api/workflows/Test/badge.svg?branch=main)

Run:

```python
pytest
```

in the root of this project.

I bootstrap the tests by getting datasets from `/app/tests/data` and from the docstrings for each path. The output is checked against the corresponding output in `/app/tests/output` and the docstring Sample Output.

I don't do any "smart" testing (e.g. test if `float`s are close to each other). Instead of writing a specific test for each endpoint I thought it was easier to use the datasets provided by Rosalind. To deal with the float issue, I round all `float`s in the outputs to 4 decimal places and then truncate to 3 decimal places. That way I can compare the rounded floats found in the docstrings to actual results. I don't want to do any rounding when I return from my endpoints.
