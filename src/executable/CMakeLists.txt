add_executable(cylaks.exe "main.cpp")
target_link_libraries(cylaks.exe cylaks)
install(TARGETS cylaks.exe DESTINATION bin)