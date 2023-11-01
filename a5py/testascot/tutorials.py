import os
import glob
import subprocess

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError

# Stores results in existing notebook
inplace = True

os.chdir("../../doc/tutorials/")
notebooks = glob.glob("*.ipynb")

# Kernel must be specified when using a virtual environment
#preprocessor = ExecutePreprocessor(timeout=600, kernel="feature8")
preprocessor = ExecutePreprocessor(timeout=600)
errors = {}
for nb in ["distributions.ipynb"]:
    #if nb == "slowingdown.ipynb": continue
    subprocess.run(["rm", "-f", "ascot.h5"])
    with open(nb) as f:
        nbin = nbformat.read(f, nbformat.NO_CONVERT)
    try:
        nbout,_ = preprocessor.preprocess(nbin)
        if inplace:
            with open(nb, "w", encoding="utf-8") as f:
                nbformat.write(nbout, f)
    except CellExecutionError as err:
        errors[nb] = err

subprocess.run(["rm", "-f", "ascot.h5"])
if len(errors) > 0:
    print("\nFollowing exceptions were logged:\n")
    for f,e in errors.items():
        print(f + ":\n")
        print(e)
        print("")
    raise Exception("Failed to run tutorials.")

print("Tutorials completed without errors")
