echo "Please select which demo you would like to run: "
echo "  1. End-tag formation"
echo "  2. Tubulin heterogeneity"
echo "  3. Kinesin heterodimer"
echo "  4. Filament separation"

read MODE 

case "$MODE" in 

    1)
        ./sim params/params_endtag.yaml demo_endtag 
        ;;
    2)
        ./sim params/params_endtag.yaml demo_tubulin hetero_tubulin
        ;;
    3)
        ./sim params/params_heterodimer.yaml demo_heterodimer kinesin_mutant
        ;;
    4)
        ./sim params/params_separation.yaml demo_separation filament_separation
        ;;

esac 