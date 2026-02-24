#!/bin/bash

# A script to run the TeraFlux Python analysis iteratively for multiple input files.
# It loops from a starting number to an ending number, generating the
# necessary filenames for each iteration.

echo "Starting TeraFlux iterative simulations..."

# Define the directory where results will be saved
RESULTS_DIR="./results:"

# Create the results directory if it doesn't already exist
# The '-p' flag ensures no error is thrown if the directory already exists.
mkdir -p "$RESULTS_DIR"

# Loop for i from 0 to 999.
# The sequence {0..999} generates all numbers in this range.
for i in {0..999}
do
  # Define the specific filenames for this iteration using the loop variable 'i'
  INPUT_FILE="inputFiles/inputData_${i}.csv"
  PREFIX_LOG_FILE="teraflux_iJO1366_${i}"

  # Add a message to the console to track progress
  echo "----------------------------------------------------"
  echo "Running iteration: $i"
  echo "Input file: $INPUT_FILE"
  echo "----------------------------------------------------"

  # Check if the specific input file for this iteration exists before running
  if [ -f "$INPUT_FILE" ]; then
    # Execute the Python script with the correct command-line arguments.
    # We assume your Python script is named 'optimized_teraflux.py'.
    # Adjust the name 'python3' if you use a different command (e.g., 'python').
    python3 ../../teraflux.py "$INPUT_FILE" -o "$RESULTS_DIR" -p "$PREFIX_LOG_FILE"
  else
    # Print a warning if an input file is missing
    echo "Warning: Input file '$INPUT_FILE' not found. Skipping iteration $i."
  fi
done

echo "All iterations complete."

