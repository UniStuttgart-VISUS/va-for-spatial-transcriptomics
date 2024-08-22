#!/bin/bash
# activate virtual environment
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    source ./venv/Scripts/activate
else
    source ./venv/bin/activate
fi
# start server
bokeh serve --show server_app.py

