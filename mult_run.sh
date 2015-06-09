#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for I in  0.000001 #0.000007 0.00001 0.00003 0.00004 0.00005 0.00006 0.00007 0.00008
  do
	 ./Release/Sibernetic -test timestep=$I timelimit=0.05
 done
