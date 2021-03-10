set(EXECUTABLE_OUTPUT_PATH ${CyLaKS_SOURCE_DIR})

file(GLOB HEADER_LIST CONFIGURE_DEPENDS "*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "*.cpp")

add_executable(cylaks ${SOURCE_LIST})
target_compile_features(cylaks PRIVATE cxx_std_17)
target_include_directories(cylaks PUBLIC ${GSL_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIRS})
target_link_libraries(cylaks yaml-cpp GSL::gsl)