## Set include directories
include_directories (BEFORE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

## Find packages
find_package(NumPy)

message (STATUS ${NUMPY_INCLUDES})
message (STATUS ${NUMPY_INCLUDES})

if (NOT PYTHONLIBS_FOUND)
    find_package(PythonLibs)
endif (NOT PYTHONLIBS_FOUND)

## Find Python site-packages directory for installation
if (NOT PYTHON_SITE_PACKAGES)
    execute_process (COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
endif (NOT PYTHON_SITE_PACKAGES)

## Determine if build will work
if (K3MATCH_WITH_PYTHON_BINDINGS AND NOT PYTHONLIBS_FOUND)
    message (WARNING "[K3MATCH] Python bindings disabled, cannot find Python libraries")
    set (K3MATCH_WITH_PYTHON_BINDINGS OFF)
endif (K3MATCH_WITH_PYTHON_BINDINGS AND NOT PYTHONLIBS_FOUND)

if (K3MATCH_WITH_PYTHON_BINDINGS AND K3MATCH_SINGLE_PRECISION)
    message (WARNING "[K3MATCH] Python bindings disabled, cannot build against single precision library.")
    set (K3MATCH_WITH_PYTHON_BINDINGS OFF)
endif (K3MATCH_WITH_PYTHON_BINDINGS AND K3MATCH_SINGLE_PRECISION)

## Add target for build and installation
if (K3MATCH_WITH_PYTHON_BINDINGS)
    message (STATUS "------------------------------------------------------------")
    message (STATUS "[PyK3MATCH] Configuration summary.")
    message (STATUS "------------------------------------------------------------")
    message (STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
    message (STATUS "PYTHON_INCLUDE_PATH = ${PYTHON_INCLUDE_PATH}")
    message (STATUS "PYTHON_SITE_PACKAGES = ${PYTHON_SITE_PACKAGES}")
    message (STATUS "NUMPY_INCLUDES = ${NUMPY_INCLUDES}")
    message (STATUS "------------------------------------------------------------")

    include_directories(${PYTHON_INCLUDE_PATH})
    include_directories(${NUMPY_INCLUDES})

    add_library (pyk3match SHARED k3match.c)

    target_link_libraries (pyk3match k3match m ${PYTHON_LIBRARIES})

    set_target_properties(pyk3match
        PROPERTIES
        PREFIX ""
        OUTPUT_NAME "k3match"
        SUFFIX ".so"
        )

    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/pyk3match.so ${CMAKE_CURRENT_BINARY_DIR}/pyk3match.py
             DESTINATION ${PYTHON_SITE_PACKAGES})
endif (K3MATCH_WITH_PYTHON_BINDINGS)

