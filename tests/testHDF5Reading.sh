#!/bin/bash
set -euo pipefail

python3 ./makeInput.py
./testHDF5Reading
