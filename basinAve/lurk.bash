#!/bin/bash

job=11781
#echo $job

sleep 3600

bash ./bat.bash
status=`squeue -u jgpender | grep $job |  tail -c 29 | cut -d " " -f1`
#echo $status

if [ $status == "R" ]
then
	bash ./lurk.bash &
fi
