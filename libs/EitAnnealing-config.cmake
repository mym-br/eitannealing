get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/lib/eitannealinglib-targets.cmake)
get_filename_component(EitAnnealing_INCLUDE_DIRS "${SELF_DIR}/include" ABSOLUTE)
