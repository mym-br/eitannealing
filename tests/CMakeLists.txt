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

add_executable(bench_symmetric_x_full bench_symmetric_x_full.cpp $<TARGET_OBJECTS:eit_time> $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET bench_symmetric_x_full PROPERTY CXX_STANDARD 17)
set_property(TARGET bench_symmetric_x_full PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_eigen_complex test_eigen_complex.cpp $<TARGET_OBJECTS:eit_time> $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET test_eigen_complex PROPERTY CXX_STANDARD 17)
set_property(TARGET test_eigen_complex PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_lb_solution test_lb_solution.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET test_lb_solution PROPERTY CXX_STANDARD 17)
set_property(TARGET test_lb_solution PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_lb_solver test_lb_solver.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
set_property(TARGET test_lb_solver PROPERTY CXX_STANDARD 17)
set_property(TARGET test_lb_solver PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(bench_new_sparsecomplexqr bench_new_sparsecomplexqr.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine> $<TARGET_OBJECTS:eit_time>)
set_property(TARGET bench_new_sparsecomplexqr PROPERTY CXX_STANDARD 17)
set_property(TARGET bench_new_sparsecomplexqr PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(bench_lbsolver_scalar_qr_vs_llt bench_lbsolver_scalar_qr_vs_llt.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine> $<TARGET_OBJECTS:eit_time>)
set_property(TARGET bench_lbsolver_scalar_qr_vs_llt PROPERTY CXX_STANDARD 17)
set_property(TARGET bench_lbsolver_scalar_qr_vs_llt PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_lb_least_eigenvalue test_lb_least_eigenvalue.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine> $<TARGET_OBJECTS:eit_time>)
set_property(TARGET test_lb_least_eigenvalue PROPERTY CXX_STANDARD 17)
set_property(TARGET test_lb_least_eigenvalue PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(test_heap test_heap.cpp $<TARGET_OBJECTS:eit_time>)
set_property(TARGET test_heap PROPERTY CXX_STANDARD 17)
set_property(TARGET test_heap PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(bench_build_qr bench_build_qr.cpp $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine> $<TARGET_OBJECTS:eit_time>)
set_property(TARGET bench_build_qr PROPERTY CXX_STANDARD 17)
set_property(TARGET bench_build_qr PROPERTY CXX_STANDARD_REQUIRED ON)
