#!/bin/bash

echo "Generating output"
./write-hdf5-output write-hdf5-output-input.ini 2>&1 | tee HDF5_writing.log

# Clean up
rm -f write-hdf5-output*.hdf5 write-hdf5-output.xmf HDF5_writing.log 

echo "Test complete"
