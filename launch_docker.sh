#! /bin/bash

# Shamelessly stolen from jeffmm/simcore project

# Run to launch a docker container named "cylaks".
# Provide the flag '-b' to force rebuilding the project image:
# >> ./run_docker.sh -b
# Provide the flag '-e' to force rebuilding the project image environment:
# >> ./run_docker.sh -e

show_help() {
    echo "Without additional options, $0 launches a docker container named 'cylaks' to run in the background"
    echo "USAGE:"
    echo "  $0 [-hbe]"
    echo "OPTIONS:"
    echo "  -h      show this menu"
    echo "  -b      force rebuild of CyLaKS docker image"
    echo "  -e      force rebuild of CyLaKS base environment docker image"
}

# Reset in case getopts has been used previously in the shell.
OPTIND=1 
build=false
build_base=false
while getopts "h?be" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
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

if $build; then
    echo "Building CyLaKS docker image"
    docker build --no-cache -t shfiorenza/cylaks:latest docker
elif $build_base; then
    echo "Building CyLaKS base environment docker image"
    docker build --no-cache -f docker/Dockerfile_base -t shfiorenza/cylaks_base:latest docker
else
    echo "Launching CyLaKS docker container"
    docker run --rm -itd -v "${PWD}":/mnt --name "cylaks" shfiorenza/cylaks bash
fi