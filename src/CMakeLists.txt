# Shamelessly stolen from jeffmm/simcore project

# find required libraries
find_package(GSL REQUIRED)
find_package(yaml-cpp REQUIRED)
message("-- Found yaml-cpp: ${YAML_CPP_INCLUDE_DIR}")

set(LIB ${GSL_LIBRARIES} ${YAML_CPP_LIBRARIES})
set(INCLUDES ${GSL_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIR})
include_directories(${INCLUDES})

file(GLOB HEADERS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/cylaks/*.hpp")
file(GLOB SOURCES "*.cpp")

# Make a static library (tradeoff of memory for better performance)
add_library(cylaks STATIC ${SOURCES} ${HEADERS})
# Add an alias so the target can be used while building
add_library(cylaks::cylaks ALIAS cylaks)
target_link_libraries(cylaks ${LIB})

# Set the include directory for builds (tests, examples, etc) and for clients
target_include_directories(cylaks PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Associate the include files with the target, so they will be installed
set_target_properties(cylaks PROPERTIES PUBLIC_HEADER "${HEADERS}")

# Set install targets and destinations
install(TARGETS cylaks
    EXPORT cylaksTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/cylaks)

# Export additional install files (for clients using CMake)
install(EXPORT cylaksTargets
    NAMESPACE cylaks::
    FILE cylaksConfig.cmake
    DESTINATION lib/cmake/cylaks)

# Add the executable subdirectory
add_subdirectory("executable")