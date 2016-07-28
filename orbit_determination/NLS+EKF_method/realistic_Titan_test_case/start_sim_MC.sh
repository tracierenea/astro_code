#!/bin/bash

# Clear the terminal screen.
clear;

# Check for a valid argument; currently there are only 3 test cases.
num_args=$#
if [ $num_args != 2 ]; then
  echo "Give two arguments: which test case to evaluate (1, 2, or 3)"
  echo "and how many Monte Carlo runs to do."
  echo -e "\nFor example:\n > bash start_sim 1 100\n"
  exit
fi

if [[ $1 != 1 ]] && [[ $1 != 2 ]]  && [[ $1 != 3 ]]; then
  echo "First argument is invalid. Choose test case 1, 2, or 3."
  echo -e "For example:\n > bash start_sim 1 100\n"
  exit
fi
echo -e "Running script for test case $1.\n"

# Run the main octave file for the specified test case.
echo "Executing command: "
echo "  > octave --silent --eval main_MC($1, $2)"
octave --silent --eval "main_MC($1, $2)"

# Create plots using python script.
python create_plots_MC.py $1;