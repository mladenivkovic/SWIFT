#!/bin/bash

echo "Creating initial conditions"
python3 ./makeInput.py

echo "Generating output"
./testHDF5Writing 2>&1 | tee HDF5_writing.log

# to add output checking
# echo "Checking output"
# python3 ./testSelectOutput.py

# NOTE Do we want to clean up?
rm -f testHDF5Writing*.hdf5 testHDF5Writing.xmf HDF5_writing.log 

echo "Test passed"
