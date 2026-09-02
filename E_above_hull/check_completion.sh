#!/bin/bash

jobs=$(squeue --me | awk '/lwalterb/ {print $1}')
for job in $jobs; do
	compl=$(cat results.out.$job | grep 'Finished job' | wc -l)
	
	mlip=$(cat results.out.$job | awk '/batchnumber 0/ {print $2}')
	error=$(grep -o 'Failed' ExpensiveMLIPs/mpid_energy_dict_$mlip.json | wc -l)

	rtime=$(squeue -j $job | awk '/lwalterb/ {print $6}')

	echo "Run time:" $rtime

	echo $job $mlip commpleted: $compl errors: $error
	echo ""
	tail -n 12 results.out.$job
	echo ""
	echo ""
done
