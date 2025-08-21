SBATCH_HEADER_FILE='scripts/sbatch_large.txt'
R1260_GROUPS=( "417" "156" "466" "272" "412" "189" "234" "485" "052" "391" 
               "294" "450" "139" "349" "462" "167" "183" "077" "481" "304" 
               "241" "338" "110" "006" "028" "121" 
               "991" "992" "993" "994" "995" "996")

align_method="all_heavy_atom"
for neighbourhood_dist in 6 10 12 20; do 
  for GROUP in "${R1260_GROUPS[@]}"; do 
    JOB_NAME="g${GROUP}_${neighbourhood_dist}_${align_method}" 
    
    # Check if a job with the same name is already running
    if squeue -n "${JOB_NAME}" --format="%j" | grep ${JOB_NAME}; then
        echo "Job ${JOB_NAME} is already running. Skipping..."
        continue
    fi
    # avg_maps/neighborhoods_${neighbourhood_dist}_${align_method}/sample5/*scatter.mrc
    # Check if there are already 32 files in the specified folder
    num_files=$(ls avg_maps/neighborhoods_${neighbourhood_dist}_${align_method}/sample?/${GROUP}*scatter.mrc avg_maps/neighborhoods_${neighbourhood_dist}_${align_method}/${GROUP}*scatter.mrc 2>/dev/null | wc -l)
    if [ "${num_files}" -ge "36" ]; then
        echo "Already ${num_files} files for GROUP=${GROUP} in ${JOB_NAME}. Skipping..."
        continue
    fi

    # Create a SBATCH job script
    {
      cat $SBATCH_HEADER_FILE
      echo "echo \"Running for GROUP=${GROUP}, neighbourhood_dist=${neighbourhood_dist}\""

      # Loop through residues and add to the job script
        echo "python scripts/3-casp16_water-get_density.py \\"
        echo "    --map scripts/Con2-2.2A_sh.mrc \\"
        echo "    --group $GROUP \\"
        echo "    --neighborhoods neighborhoods_${neighbourhood_dist} \\"
        echo "    --output_folder avg_maps  \\"
        echo "    --alignment_method $align_method \\"
    } > "setup_sbatch/${JOB_NAME}.sh"
    sed -i "s/XXX/${JOB_NAME}/g" "setup_sbatch/${JOB_NAME}.sh"
    sbatch "setup_sbatch/${JOB_NAME}.sh"

    for ((i=1; i<=5; i++)); do
      JOB_NAME_EXT="${JOB_NAME}_s${i}"
      # Check if a job with the same name is already running
      if squeue -n "${JOB_NAME_EXT}" --format="%j" | grep ${JOB_NAME_EXT}; then
        echo "Job ${JOB_NAME_EXT} is already running. Skipping..."
        continue
      fi

      # Create a SBATCH job script
      {
        cat $SBATCH_HEADER_FILE
        echo "echo \"Running for GROUP=${GROUP}, neighbourhood_dist=${neighbourhood_dist}, sample=${i}\""

        # Add the python command to the job script
        echo "python scripts/3-casp16_water-get_density.py \\"
        echo "    --map scripts/Con2-2.2A_sh.mrc \\"
        echo "    --group $GROUP \\"
        echo "    --neighborhoods neighborhoods_${neighbourhood_dist} \\"
        echo "    --output_folder avg_maps  \\"
        echo "    --alignment_method $align_method \\"
        echo "    --sample $i \\"
      } > "setup_sbatch/${JOB_NAME_EXT}.sh"

      sed -i "s/XXX/${JOB_NAME_EXT}/g" "setup_sbatch/${JOB_NAME_EXT}.sh"
      sbatch "setup_sbatch/${JOB_NAME_EXT}.sh"
    done
  done
done

neighbourhood_dist=10
align_method="all_heavy_atom"
for align_method in "3_atom" "5_atom" "backbone"; do 
  for GROUP in "${R1260_GROUPS[@]}"; do 
    JOB_NAME="g${GROUP}_${neighbourhood_dist}_${align_method}" 
    
    # Check if a job with the same name is already running
    if squeue -n "${JOB_NAME}" --format="%j" | grep ${JOB_NAME}; then
        echo "Job ${JOB_NAME} is already running. Skipping..."
        continue
    fi
    # avg_maps/neighborhoods_${neighbourhood_dist}_${align_method}/sample5/*scatter.mrc
    # Check if there are already 32 files in the specified folder
    num_files=$(ls avg_maps/neighborhoods_${neighbourhood_dist}_${align_method}/sample?/${GROUP}*scatter.mrc avg_maps/neighborhoods_${neighbourhood_dist}_${align_method}/${GROUP}*scatter.mrc 2>/dev/null | wc -l)
    if [ "${num_files}" -ge "36" ]; then
        echo "Already ${num_files} files for GROUP=${GROUP} in ${JOB_NAME}. Skipping..."
        continue
    fi

    # Create a SBATCH job script
    {
      cat $SBATCH_HEADER_FILE
      echo "echo \"Running for GROUP=${GROUP}, neighbourhood_dist=${neighbourhood_dist}\""

      # Loop through residues and add to the job script
        echo "python scripts/3-casp16_water-get_density.py \\"
        echo "    --map scripts/Con2-2.2A_sh.mrc \\"
        echo "    --group $GROUP \\"
        echo "    --neighborhoods neighborhoods_${neighbourhood_dist} \\"
        echo "    --output_folder avg_maps  \\"
        echo "    --alignment_method $align_method \\"
    } > "setup_sbatch/${JOB_NAME}.sh"
    sed -i "s/XXX/${JOB_NAME}/g" "setup_sbatch/${JOB_NAME}.sh"
    sbatch "setup_sbatch/${JOB_NAME}.sh"

    for ((i=1; i<=5; i++)); do
      JOB_NAME_EXT="${JOB_NAME}_s${i}"
      # Check if a job with the same name is already running
      if squeue -n "${JOB_NAME_EXT}" --format="%j" | grep ${JOB_NAME_EXT}; then
        echo "Job ${JOB_NAME_EXT} is already running. Skipping..."
        continue
      fi

      # Create a SBATCH job script
      {
        cat $SBATCH_HEADER_FILE
        echo "echo \"Running for GROUP=${GROUP}, neighbourhood_dist=${neighbourhood_dist}, sample=${i}\""

        # Add the python command to the job script
        echo "python scripts/3-casp16_water-get_density.py \\"
        echo "    --map scripts/Con2-2.2A_sh.mrc \\"
        echo "    --group $GROUP \\"
        echo "    --neighborhoods neighborhoods_${neighbourhood_dist} \\"
        echo "    --output_folder avg_maps  \\"
        echo "    --alignment_method $align_method \\"
        echo "    --sample $i \\"
      } > "setup_sbatch/${JOB_NAME_EXT}.sh"

      sed -i "s/XXX/${JOB_NAME_EXT}/g" "setup_sbatch/${JOB_NAME_EXT}.sh"
      sbatch "setup_sbatch/${JOB_NAME_EXT}.sh"
    done
  done
done
