 ###############################################################################
 # This file is part of SWIFT.
 # Copyright (C) 2026 Elijah Cann (e.a.cann72@gmail.com).
 #
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 #
 ###############################################################################

noinst_LTLIBRARIES = libwritehdf5output.la

libwritehdf5output_la_CFLAGS = -I$(top_srcdir)/src $(HDF5_CPPFLAGS)

libwritehdf5output_la_SOURCES = \
    $(top_srcdir)/src/threadpool.c \
    $(top_srcdir)/src/clocks.c \
    $(top_srcdir)/src/parser.c \
    $(top_srcdir)/src/output_options.c \
    $(top_srcdir)/src/units.c \
    $(top_srcdir)/src/physical_constants.c \
    $(top_srcdir)/src/single_io.c \
    $(top_srcdir)/src/cosmology.c \
    $(top_srcdir)/src/hydro_properties.c \
    $(top_srcdir)/src/common_io_fields.c \
    $(top_srcdir)/src/lightcone/lightcone_map_types.c \
    $(top_srcdir)/src/common_io.c \
    $(top_srcdir)/src/common_io_copy.c \
    $(top_srcdir)/src/ic_info.c \
    $(top_srcdir)/src/output_list.c \
    $(top_srcdir)/src/gravity_properties.c \
    $(top_srcdir)/src/hdf5_object_to_blob.c \
    $(top_srcdir)/src/tools.c \
    $(top_srcdir)/src/version.c \
    $(top_srcdir)/src/restart.c \
    $(top_srcdir)/src/part_type.c \
    $(top_srcdir)/src/io_compression.c \
    $(top_srcdir)/src/xmf.c \
    $(top_srcdir)/src/common_io_cells.c \
    $(top_srcdir)/src/equation_of_state.c \
    write-hdf5-output-stubs.c


noinst_PROGRAMS = write-hdf5-output
write_hdf5_output_SOURCES = write-hdf5-output.c
write_hdf5_output_CFLAGS  = -I$(top_srcdir)/src $(HDF5_CPPFLAGS) $(GSL_INCS)
write_hdf5_output_LDADD   = libwritehdf5output.la
write_hdf5_output_LDFLAGS = $(HDF5_LDFLAGS) $(HDF5_LIBS) $(GSL_LIBS) -lm

EXTRA_DIST = write-hdf5-output.sh HDF5WritingParameters.yml HDF5input.hdf5
