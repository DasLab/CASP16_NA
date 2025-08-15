# get the neighborhoods
# note this was done on a cpu cluster, SBATCH_HEADER_FILE
# is just a sbatch header file for running each job
SBATCH_HEADER_FILE='scripts/sbatch.txt'
align_method="all_heavy_atom"
for neighbourhood_dist in 6 10 12 20; do 
  for RES in {22..408}; do 
    JOB_NAME="${RES}_${neighbourhood_dist}_${align_method}" 
    
    # Check if a job with the same name is already running
    if squeue -n "${JOB_NAME}" --format="%j" | grep ${JOB_NAME}; then
        echo "Job ${JOB_NAME} is already running. Skipping..."
        continue
    fi

    # Check if there are already 32 files in the specified folder
    num_files=$(ls neighborhoods_${neighbourhood_dist}/${align_method}/*/${RES}_rna.pdb 2>/dev/null | wc -l)
    if [ "${num_files}" -ge "32" ]; then
        echo "Already ${num_files} files for RES=${RES} in ${neighborhood_folder}. Skipping..."
        continue
    fi

    # Create a SBATCH job script
    {
      cat $SBATCH_HEADER_FILE 
      echo "echo \"Running for RES=${RES}, neighbourhood_dist=${neighbourhood_dist}\""

      # Loop through residues and add to the job script
        echo "python scripts/2-casp16_water-get_areas_of_interest.py \\"
        echo "    --native_csv scripts/9CBU-rna.csv \\"
        echo "    --input_folder R1260_serial \\"
        echo "    --neighborhood_output neighborhoods_${neighbourhood_dist} \\"
        echo "    --res $RES  \\"
        echo "    --neighborhood_radius ${neighbourhood_dist} \\"
        echo "    --solvent_radius 5 \\"
        echo "    --alignment_method $align_method \\"
        echo "    --num_samples 5 "
    } > "setup_sbatch/${JOB_NAME}.sh"

    sed -i "s/XXX/${JOB_NAME}/g" "setup_sbatch/${JOB_NAME}.sh"
    sbatch "setup_sbatch/${JOB_NAME}.sh"
  done
done

neighbourhood_dist=10
for align_method in "3_atom" "5_Atom" "backbone"; do 
  for RES in {22..408}; do 
    JOB_NAME="${RES}_${neighbourhood_dist}_${align_method}" 
    
    # Check if a job with the same name is already running
    if squeue -n "${JOB_NAME}" --format="%j" | grep ${JOB_NAME}; then
        echo "Job ${JOB_NAME} is already running. Skipping..."
        continue
    fi

    # Check if there are already 32 files in the specified folder
    num_files=$(ls neighborhoods_${neighbourhood_dist}/${align_method}/*/${RES}_rna.pdb 2>/dev/null | wc -l)
    if [ "${num_files}" -ge "32" ]; then
        echo "Already ${num_files} files for RES=${RES} in ${neighborhood_folder}. Skipping..."
        continue
    fi

    # Create a SBATCH job script
    {
      cat $SBATCH_HEADER_FILE 
      echo "echo \"Running for RES=${RES}, neighbourhood_dist=${neighbourhood_dist}\""

      # Loop through residues and add to the job script
        echo "python scripts/2-casp16_water-get_areas_of_interest.py \\"
        echo "    --native_csv scripts/9CBU-rna.csv \\"
        echo "    --input_folder R1260_serial \\"
        echo "    --neighborhood_output neighborhoods_${neighbourhood_dist} \\"
        echo "    --res $RES  \\"
        echo "    --neighborhood_radius ${neighbourhood_dist} \\"
        echo "    --solvent_radius 5 \\"
        echo "    --alignment_method $align_method \\"
        echo "    --num_samples 5 "
    } > "setup_sbatch/${JOB_NAME}.sh"

    sed -i "s/XXX/${JOB_NAME}/g" "setup_sbatch/${JOB_NAME}.sh"
    sbatch "setup_sbatch/${JOB_NAME}.sh"
  done
done
