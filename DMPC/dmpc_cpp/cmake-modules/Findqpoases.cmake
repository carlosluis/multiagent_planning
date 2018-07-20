# Find qpOASES header and library.
#

# This module defines the following uncached variables:
#  QPOASES_FOUND, if false, do not try to use qpOASES.
#  QPOASES_INCLUDE_DIRS, where to find qpOASES/include
#  QPOASES_LIBRARIES, the libraries to link against to use the qpoases library
#  QPOASES_LIBRARY_DIRS, the directory where the qpoases library is found.

# look for QPOASES
find_library(QPOASES_LIBRARY
        NAMES qpOASES
        HINTS $ENV{HOME}/qpOASES-3.2.1/bin)
message("QPoases library: ${QPOASES_LIBRARY}")
find_path(QPOASES_INCLUDE_DIR qpOASES PATHS $ENV{HOME}/qpOASES-3.2.1/include)
message("QPoases library: ${QPOASES_INCLUDE_DIR}")


set(QPOASES_LIBRARY_DIR "")
get_filename_component(QPOASES_LIBRARY_DIRS ${QPOASES_LIBRARY} PATH)
# Set uncached variables as per standard.
set(QPOASES_FOUND ON)
set(QPOASES_INCLUDE_DIRS ${QPOASES_INCLUDE_DIR})
set(QPOASES_LIBRARIES ${QPOASES_LIBRARY})


if(QPOASES_LIBRARY)
  if(NOT QPOASES_FIND_QUIETLY)
    message(STATUS "FindQPOASES: Found both include/qpoases and libqpOASES.a")
  endif(NOT QPOASES_FIND_QUIETLY)
else(QPOASES_LIBRARY)
  if(QPOASES_INCLUDE_DIR)
    message(FATAL_ERROR "FindQPOASES: Could not find libqpoases.a")
  endif(QPOASES_INCLUDE_DIR)
endif(QPOASES_LIBRARY)

