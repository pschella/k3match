## Set include directories
include_directories (BEFORE ${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

## List of source files
FILE (GLOB k3match_sources *.c)

## List of header files
FILE (GLOB k3match_headers *.h)

## Build the library
if (K3MATCH_SHARED_LIBRARY)
    add_library (k3match SHARED ${k3match_sources})
else (K3MATCH_SHARED_LIBRARY)
    add_library (k3match STATIC ${k3match_sources})
endif (K3MATCH_SHARED_LIBRARY)

## Library versioning
set_target_properties (k3match PROPERTIES VERSION ${K3MATCH_VERSION})

## Install libray
install (TARGETS k3match
  DESTINATION lib
  )

## Install header files
install (FILES ${k3match_headers}
  DESTINATION include/k3match
  )

