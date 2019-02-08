add_executable(test_incomplete_qr test_incomplete_lq.cpp $<TARGET_OBJECTS:eit_core>)
set_property(TARGET test_incomplete_qr PROPERTY CXX_STANDARD 17)
set_property(TARGET test_incomplete_qr PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(stddevtest stddevtest.cpp)
set_property(TARGET stddevtest PROPERTY CXX_STANDARD 17)
set_property(TARGET stddevtest PROPERTY CXX_STANDARD_REQUIRED ON)

add_test(incomplete_qr COMMAND test_incomplete_qr)
add_test(stdde COMMAND stddevtest)
