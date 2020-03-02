add_executable(test_incomplete_qr test_incomplete_lq.cpp $<TARGET_OBJECTS:eit_core>)
set_property(TARGET test_incomplete_qr PROPERTY CXX_STANDARD 17)
set_property(TARGET test_incomplete_qr PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(stddevtest stddevtest.cpp)
set_property(TARGET stddevtest PROPERTY CXX_STANDARD 17)
set_property(TARGET stddevtest PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_prereg_qr test_prereg_complex_incomplete.cpp $<TARGET_OBJECTS:eit_core>)
set_property(TARGET test_prereg_qr PROPERTY CXX_STANDARD 17)
set_property(TARGET test_prereg_qr PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_solver_lb_complex test_solver_lb_complex.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET test_solver_lb_complex PROPERTY CXX_STANDARD 17)
set_property(TARGET test_solver_lb_complex PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(bench_symmetric_x_full bench_symmetric_x_full.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET bench_symmetric_x_full PROPERTY CXX_STANDARD 17)
set_property(TARGET bench_symmetric_x_full PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_eigen_complex test_eigen_complex.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET test_eigen_complex PROPERTY CXX_STANDARD 17)
set_property(TARGET test_eigen_complex PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_lb_solution test_lb_solution.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET test_lb_solution PROPERTY CXX_STANDARD 17)
set_property(TARGET test_lb_solution PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_lb_solver test_lb_solver.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET test_lb_solver PROPERTY CXX_STANDARD 17)
set_property(TARGET test_lb_solver PROPERTY CXX_STANDARD_REQUIRED ON)
