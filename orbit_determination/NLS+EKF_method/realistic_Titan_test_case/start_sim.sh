#!/bin/bash

# Clear the terminal screen.
clear;

# Check for a valid argument; currently there are only 3 test cases.
num_args=$#
if [ $num_args != 1 ]; then
  echo "Give one argument for which test case to evaluate: 1,2, or 3"
  echo -e "For example:\n > bash start_sim 1\n"
  exit
fi

if [[ $1 != 1 ]] && [[ $1 != 2 ]]  && [[ $1 != 3 ]]; then
  echo "Invalid argument. Choose test case 1, 2, or 3."
  echo -e "For example:\n > bash start_sim 1\n"
  exit
fi
echo -e "Running script for test case $1.\n"

# These are the files containing simulation data.
file1=measurement_data.txt
file2=sat_states.txt

# Remove data files from the previous run.
echo "Checking for data files from a previous run:"
if [ -f "$file1" ];
then
  echo "  - Deleting file $file1."
  rm $file1
else
  echo "  - File $file1 does not exist"
fi

if [ -f "$file2" ];
then
  echo "  - Deleting file $file2."
  rm $file2
else
  echo "  - File $file2 does not exist"
fi
echo ""

# Run the main octave file for the specified test case.
echo "Executing command: "
echo "  > octave --silent --eval main($1)"
octave --silent --eval "main($1)"

# Create plots using python script.
python create_plots.py  $1;
# if [ -f "$file1" ] && [ -f "$file2" ];
# then
#   python create_plots.py  $1;
# else
#   echo -e "\nNo data. Skipping plots.\n"
# fi