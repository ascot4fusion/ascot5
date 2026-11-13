"""
Execute all tutorials and log any errors.

Exception is raised in case any of the tutorials fail. Otherwise results are
written in place and this script exits gracefully.
"""

import glob

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError

store_results_in_existing_notebook = True
notebooks = glob.glob("*.ipynb")
preprocessor = ExecutePreprocessor(timeout=600)
errors = {}
for nb in notebooks:
    print(f"Processing notebook {nb}")
    with open(nb) as f:
        nbin = nbformat.read(f, nbformat.NO_CONVERT)
    try:
        nbout, _ = preprocessor.preprocess(nbin)
        if store_results_in_existing_notebook:
            with open(nb, "w", encoding="utf-8") as f:
                nbformat.write(nbout, f)
    except CellExecutionError as err:
        errors[nb] = err

if len(errors) > 0:
    print("\nFollowing exceptions were logged:\n")
    for f, e in errors.items():
        print(f + ":\n")
        print(e)
        print("")
    raise Exception("One or more tutorials .")

print("Tutorials completed.")
