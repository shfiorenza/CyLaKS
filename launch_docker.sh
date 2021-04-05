#! /bin/bash

# Shamelessly stolen from jeffmm/simcore project

# Run to launch a docker container named "cylaks".
# Provide the flag '-b' to force rebuilding the project image:
# >> ./run_docker.sh -b
# Provide the flag '-e' to force rebuilding the project image environment:
# >> ./run_docker.sh -e
# Provide the flag '-d' to to launch run_demos.sh script in Docker
# >> ./run_docker.sh -d

show_help() {
    echo "Without additional options, $0 launches a docker container named 'cylaks' to run in the background"
    echo "USAGE:"
    echo "  $0 [-hbed]"
    echo "OPTIONS:"
    echo "  -h      show this menu"
    echo "  -b      force rebuild of CyLaKS docker image"
    echo "  -e      force rebuild of CyLaKS base environment docker image"
    echo "  -d      launch run_demos.sh script in Docker container"
}

# Reset in case getopts has been used previously in the shell.
OPTIND=1 
build=false
build_base=false
launch_demos=false
while getopts "h?bed" opt; do
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
    d)
        launch_demos=true
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
elif $launch_demos; then
    echo "Please select which demo you would like to run: "
    echo "  1. End-tag formation"
    echo "  2. Tubulin heterogeneity"
    echo "  3. Kinesin heterodimer"
    echo "  4. Filament separation"
    read MODE 
    case "$MODE" in 

        1)
            docker exec -it cylaks cylaks.exe params/params_endtag.yaml demo_endtag 
            ;;
        2)
            docker exec -it cylaks cylaks.exe params/params_endtag.yaml demo_tubulin hetero_tubulin
            ;;
        3)
            docker exec -it cylaks cylaks.exe params/params_heterodimer.yaml demo_heterodimer kinesin_mutant
            ;;
        4)
            docker exec -it cylaks cylaks.exe params/params_separation.yaml demo_separation filament_separation
            ;;

    esac 
else
    echo "Launching CyLaKS docker container"
    docker run --rm -itd -v "${PWD}":/mnt --name "cylaks" shfiorenza/cylaks bash
fi