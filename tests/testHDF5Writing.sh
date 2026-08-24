#!/bin/bash

echo "Generating output"
./testHDF5Writing testHDF5Writing.ini 2>&1 | tee HDF5_writing.log

# Clean up
rm -f testHDF5Writing*.hdf5 testHDF5Writing.xmf HDF5_writing.log 

echo "Test passed"
