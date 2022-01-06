#!/bin/bash
###
# https://flask.palletsprojects.com/en/2.0.x/installation/
# python3 -m venv venv
# . venv/bin/activate
# pip install Flask
###
# https://flask.palletsprojects.com/en/2.0.x/tutorial/blog/
###
# VirtualEnv should be activated:
# . venv/bin/activate
###
export FLASK_APP=depict
export FLASK_ENV=development
#flask init-db
flask run
###
# VirtualEnv deactivate:
# deactivate
