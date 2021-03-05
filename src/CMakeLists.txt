set(EXECUTABLE_OUTPUT_PATH ${CyLaKS_SOURCE_DIR})

file(GLOB HEADER_LIST CONFIGURE_DEPENDS "*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "*.cpp")

add_executable(sim ${SOURCE_LIST})
target_compile_features(sim PRIVATE cxx_std_17)
target_include_directories(sim PUBLIC ${GSL_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIRS})
target_link_libraries(sim yaml-cpp GSL::gsl)