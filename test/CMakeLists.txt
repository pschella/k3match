set (tests test_k3match test_median test_insert)

foreach (test ${tests})

    add_executable (${test} ${test}.c)

    target_link_libraries (${test} m k3match)

endforeach(test ${tests})

