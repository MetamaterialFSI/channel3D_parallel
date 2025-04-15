#!/bin/bash

# Define the number of experiments for mode 0
num_exp=10  # Modify this number as needed

# Define the custom list of experiments for mode 3 and mode 4 (hardcoded list)
custom_list=(5 7 10 12 15)  # Modify this list as needed

# Check if the required number of arguments are passed
if [ $# -lt 1 ]; then
  echo "Usage: $0 <mode> [additional arguments...]"
  exit 1
fi

# Capture the mode and other arguments
mode=$1

# Define actions based on the mode
case $mode in
  0) # Mode 0: Compile and launch

    # prepare libraries
    module load fftw
    module load intel-oneapi-compilers
    module load openmpi

    # Navigate to the code directory and compile
    cmake -B build
    cmake --build build

    # Loop through the experiments
    for i in $(seq 1 ${num_exp}); do
      echo "Submit ctrl_exp file ${i}..."
      cd ./${i} || { echo "Failed to enter ./${i}"; exit 1; }
      sbatch launch_job.sh
      cd ../
    done
    ;;
  
  1) # Mode 1: Delete jobs using slurm commands
    if [ $# -lt 3 ]; then
      echo "Usage for mode 1: $0 1 <id_start> <id_end>"
      exit 1
    fi

    id_start=$2
    id_end=$3

    # Loop through job IDs and cancel them
    for id in $(seq ${id_start} ${id_end}); do
      echo "Cancelling job ID: $id"
      scancel $id
    done
    ;;
  
  2) # Mode 2: Resubmit jobs
    if [ $# -lt 3 ]; then
      echo "Usage for mode 2: $0 2 <exp_start> <exp_end>"
      exit 1
    fi

    exp_start=$2
    exp_end=$3

    # Loop through the experiments and resubmit them
    for i in $(seq ${exp_start} ${exp_end}); do
      echo "Resubmitting experiment $i..."
      cd ./${i} || { echo "Failed to enter ./${i}"; exit 1; }
      sbatch launch_job
      cd ../
    done
    ;;
  
  3) # Mode 3: Submit jobs based on a custom list of experiment numbers
    echo "Submitting jobs from custom list: ${custom_list[*]}"
    
    for i in "${custom_list[@]}"; do
      echo "Submit ctrl_exp file ${i}..."
      cd ./${i} || { echo "Failed to enter ./${i}"; exit 1; }
      sbatch launch_job
      cd ../
    done
    ;;
  
  4) # Mode 4: Delete jobs based on the custom list of job IDs
    echo "Cancelling jobs from custom list: ${custom_list[*]}"
    
    for id in "${custom_list[@]}"; do
      echo "Cancelling job ID: $id"
      scancel $id
    done
    ;;

  *) # Invalid mode
    echo "Invalid mode. Use 0 for compile and launch, 1 for delete, 2 for resubmit, 3 for custom list submit, or 4 for custom list delete."
    exit 1
    ;;
esac
