add_executable (test_k3match test_k3match.c)

target_link_libraries (test_k3match m k3match)

add_executable (test_median test_median.c)

target_link_libraries (test_median m k3match)

