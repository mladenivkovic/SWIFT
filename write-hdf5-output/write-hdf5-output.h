#ifndef SWIFT_WRITE_HDF5_OUTPUT_H
#define SWIFT_WRITE_HDF5_OUTPUT_H

/**
 * @brief Generate a uniform-lattice gas IC file and read it in
 *        before writing it as a snapshot via write_output_single()
 *
 * @param numberOfParticles Approximate particle count (rounded up to the
 *        nearest cube, to allow an L x L x L latice).
 * @param param_filename Path to a SWIFT-format YAML parameter file
 *        supplying the unit system, physical constants, and output
 *        selection options (see HDF5WritingParameters.yml for an example).
 *
 * @return 0 on success, non-zero on failure.
 */
int write_hdf5_output_run(int numberOfParticles, const char *param_filename);

#endif /* SWIFT_WRITE_HDF5_OUTPUT_H */
