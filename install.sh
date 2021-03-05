#!/bin/bash

# Shamelessly stolen from jeffmm/simcore project

do_build() {
    mkdir build
    cd build || exit 1
    cmake ${CMAKE_FLAGS} ..
    make -j8
    if $build_docs; then
        make docs
    fi
    if $run_tests; then
        make test
    fi
    cd ..
}

clean_cmake_files() {
    rm -rf build/CMake*
    rm -rf build/Makefile
    rm -rf build/cmake_install.cmake
    rm -rf build/src
    rm -rf build/Doxyfile
    rm -rf build/tests
}

show_help() {
    echo "USAGE:"
    echo "  $0 [-hcdtD]"
    echo "OPTIONS:"
    echo "  -h      show this menu"
    echo "  -c      clean build directory"
    echo "  -d      (INACTIVE) build Doxygen documentation"
    echo "  -t      (INACTIVE) build and run simcore unit tests"
    echo "  -D      build simcore in Debug mode"
}

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

CMAKE_FLAGS="-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
build_docs=false
run_tests=false
while getopts "h?cdtD" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    c)  
        clean_cmake_files
        exit 0
        ;;
    d)  build_docs=true
        ;;
    t)  
        run_tests=true
        CMAKE_FLAGS="${CMAKE_FLAGS} -DTESTS=TRUE"
        ;;
    D)  CMAKE_FLAGS="${CMAKE_FLAGS} -DDEBUG=TRUE"
        ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

do_build