#! /bin/bash

# *** Shamelessly stolen/adapted from jeffmm/simcore project ***

show_help() {
    echo "USAGE:"
    echo "  $0 [flag]"
    echo "OPTIONS:"
    echo "  -h      help: show this menu (default if no options are provided)"
    echo "  -i      initialize: start a new docker container named 'cylaks' using docker image" 
    echo "  -r      run: launch CyLaKS program in 'cylaks' docker container and link output to local dir"
    echo "  -c      clean up: stop 'cylaks' docker container and remove docker image from local system"
    echo "  -b      build: force local rebuild of CyLaKS docker image (use if pulled image does not work)"
    echo "  -e      build: force local rebuild of base Linux environment docker image"
}

if [[ $# -gt 1 ]]; then
    echo "Error; please select just one flag at a time"
    exit 0
fi

OPTIND=1
init=false
run=false
clean=false 
build=false
build_base=false
while getopts "h?ircbd" opt; do
    if $init OR $run OR $clean OR $build OR $build_base; then
        echo "Error; please select just one flag at a time"
        exit 0
    fi
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    i)  
        init=true
        ;;
    r)  
        run=true
        ;;
    c)  
        clean=true
        ;;
    b)  
        build=true
        ;;
    e)
        build_base=true
        ;;
    esac
done

shift $((OPTIND-1))
[ "${1:-}" = "--" ] && shift

if $init; then
    echo "Launching CyLaKS docker container named 'cylaks'"
    docker run --rm -itd -v "${PWD}":/mnt --name "cylaks" shfiorenza/cylaks bash
elif $run; then
    echo "Running CyLaKS program in docker container. Data files will be output to local directory"
    echo "(Interactive launcher only; quick-launch syntax currently not supported w/ Docker)"
    echo ""
    docker exec -it cylaks cylaks.exe
elif $clean; then 
    echo "Stopping CyLaKS container and removing docker image from local system."
    docker container stop cylaks
    docker image rm shfiorenza/cylaks
elif $build; then
    echo "Building CyLaKS docker image"
    docker build --no-cache -t shfiorenza/cylaks:latest docker
elif $build_base; then
    echo "Building CyLaKS base environment docker image"
    docker build --no-cache -f docker/Dockerfile_base -t shfiorenza/cylaks_base:latest docker
else
    show_help
fi