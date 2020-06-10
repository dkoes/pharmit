# - Find ANN
# Find the native ANN includes and library
# This module defines
#  ANN_INCLUDE_DIR, where to find ANN/ANN.h, etc.
#  ANN_LIBRARIES, the libraries needed to use ANN.
#  ANN_FOUND, If false, do not try to use ANN.
# also defined, but not for general use are
#  ANN_LIBRARY, where to find the ANN library.

include(FindPackageHandleStandardArgs)

find_path(ANN_INCLUDE_DIR ANN/ANN.h
          HINTS ${ANN_ROOT_DIR}/include
         )
find_library(ANN_LIBRARY ann
             HINTS ${ANN_ROOT_DIR}/lib /usr/lib
)

find_package_handle_standard_args(ANN DEFAULT_MSG
    ANN_LIBRARY ANN_INCLUDE_DIR)

if(ANN_FOUND)
    set(ANN_INCLUDE_DIRS ${ANN_INCLUDE_DIR})
    set(ANN_LIBRARIES ${ANN_LIBRARY})
endif()

add_library(ann STATIC IMPORTED)
#specify what and where... 
#target_include_directories(ann INTERFACE ${ANN_INCLUDE_DIR})
