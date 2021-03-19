# Shamelessly stolen from jeffmm/simcore project

set(EXECUTABLE_OUTPUT_PATH ${cylaks_SOURCE_DIR})
add_executable(cylaks.exe "main.cpp")
target_link_libraries(cylaks.exe cylaks)
install(TARGETS cylaks.exe DESTINATION bin)