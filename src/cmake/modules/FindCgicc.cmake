# Copyright (C) 2015-2016  Thomas Pircher <tehpeh-web@tty1.net>
#
# woale -- a wiki engine.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



# FindCgicc
# ---------
#
# Find libcgicc include dirs and libraries.
#
# Use this module by invoking find_package with the form:
#
#   find_package(Cgicc
#     [REQUIRED]             # Fail with error if the library is not found
#     )
#
# This module finds headers and requested component libraries.  Results are
# reported in variables:
#
#   Cgicc_FOUND             - True if headers and requested libraries were found
#   Cgicc_INCLUDE_DIRS      - include directories
#   Cgicc_LIBRARIES         - libcgicc component libraries to be linked
#
#
# The following import target is created
#
#   Cgicc::Cgicc

# Use pkg-config (if available) to find the directories of the componants.
find_package(PkgConfig QUIET)
pkg_check_modules(PKG_CGICC QUIET "libcgicc")
set(CGICC_DEFINITIONS ${PKG_CGICC_CFLAGS_OTHER})

find_path(CGICC_INCLUDE_DIR
        NAMES cgicc/Cgicc.h
        #HINTS ${PKG_CGICC_INCLUDEDIR} ${PKG_CGICC_INCLUDE_DIRS}
        HINTS ${Cgicc_ROOT_DIR}/include)
        #PATH_SUFFIXES cgicc)

find_library(CGICC_LIBRARY
        NAMES cgicc
        #HINTS ${PKG_CGICC_LIBDIR} ${PKG_CGICC_LIBRARY_DIRS})
        HINTS ${Cgicc_ROOT_DIR}/lib)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cgicc DEFAULT_MSG
        CGICC_LIBRARY CGICC_INCLUDE_DIR)

mark_as_advanced(Cgicc_FOUND CGICC_INCLUDE_DIR CGICC_LIBRARY)


if(Cgicc_FOUND AND NOT TARGET Cgicc::Cgicc)
    add_library(Cgicc::Cgicc UNKNOWN IMPORTED)
    set_property(TARGET Cgicc::Cgicc PROPERTY IMPORTED_LOCATION ${CGICC_LIBRARY})
    set_property(TARGET Cgicc::Cgicc PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CGICC_INCLUDE_DIR}")
endif()
