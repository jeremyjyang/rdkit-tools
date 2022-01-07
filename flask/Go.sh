#!/bin/bash
###
# https://flask.palletsprojects.com/en/2.0.x/tutorial/blog/
###
# We need to use Conda, not VirtualEnv, for RDKit support.
# https://softwarejargon.com/dockerizing-python-flask-app-and-conda-environment/
###
# conda create -c conda-forge -n rdktools rdkit
# conda activate rdktools
# pip install Flask
# Or:
# conda env create -f environment.yml
# conda activate rdktools
###
source $(dirname $(which conda))/../bin/activate rdktools
export FLASK_APP=depict
export FLASK_ENV=development
flask run
conda deactivate
###
