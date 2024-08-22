#!/bin/bash
# activate virtual environment
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    source ./venv/Scripts/activate
else
    source ./venv/bin/activate
fi
# install dependencies
pip install -r requirements.txt

