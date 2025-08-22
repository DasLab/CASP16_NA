
MAP_FOLDER=/data/rkretsch/CASP-assessment/water/TEST_MAPS/avg_maps/neighborhoods_10_all_heavy_atom/
R1260_GROUPS=( "417" "156" "466" "272" "412" "189" "234" "485" "052" "391" 
               "294" "450" "139" "349" "462" "167" "183" "077" "481" "304" 
               "241" "338" "110" "006" "028" "121" 
               "991" "992" "993" "994" "995" "996" )
SIGMA=10

for GROUP in ${R1260_GROUPS[@]}; do
    # chimerax --script "image_groups.cxc $MAP_FOLDER $GROUP scatter $SIGMA"
    timeout 2m chimerax --script "image_groups.cxc $MAP_FOLDER $GROUP density $SIGMA"
done
