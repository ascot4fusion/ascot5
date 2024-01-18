#!/bin/bash
# Clean output and metadata from all tutorials (when executed in this folder).
jupyter-nbconvert --ClearOutputPreprocessor.enabled=True \
		  --ClearMetadataPreprocessor.enabled=True \
		  --ClearMetadataPreprocessor.preserve_cell_metadata_mask='[("nbsphinx-thumbnail")]' \
		                    --inplace *.ipynb
