if (NOT NUMPY_FOUND)

  ## Initialize variables
  unset (NUMPY_FOUND)
  unset (NUMPY_INCLUDES)

  ## Check for installation of Python
  if (NOT PYTHON_FOUND)
    set (PYTHON_FIND_QUIETLY ${NUMPY_FIND_QUIETLY})
    find_package (PythonInterp)
  endif (NOT PYTHON_FOUND)

  ## Search for header files
  if (PYTHON_EXECUTABLE)
    ## Use Python to determine the include directory
    execute_process (
      COMMAND ${PYTHON_EXECUTABLE} -c import\ numpy\;\ print\ numpy.get_include\(\)\;
      ERROR_VARIABLE NUMPY_FIND_ERROR
      RESULT_VARIABLE NUMPY_FIND_RESULT
      OUTPUT_VARIABLE NUMPY_FIND_OUTPUT
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    ## process the output from the execution of the command
    if (NOT NUMPY_FIND_RESULT)
      set (NUMPY_INCLUDES ${NUMPY_FIND_OUTPUT})
    endif (NOT NUMPY_FIND_RESULT)
  endif (PYTHON_EXECUTABLE)
  
  ## Actions taken after completing the search
  if (NUMPY_INCLUDES)
    set (NUMPY_FOUND TRUE)
    if (NOT NUMPY_FIND_QUIETLY)
      message (STATUS "[FindNumPy] Found components for NumPy")
      message (STATUS "NUMPY_INCLUDES  = ${NUMPY_INCLUDES}")
    endif (NOT NUMPY_FIND_QUIETLY)
  else (NUMPY_INCLUDES)
    set (NUMPY_FOUND FALSE)
    if (NOT NUMPY_FIND_QUIETLY)
      message (FATAL_ERROR "[FindNumPy] Attempt to import NumPy failed!")
    endif (NOT NUMPY_FIND_QUIETLY)
  endif (NUMPY_INCLUDES)

  ## Mark advanced variables
  mark_as_advanced (
    NUMPY_INCLUDES
    )

endif (NOT NUMPY_FOUND)
