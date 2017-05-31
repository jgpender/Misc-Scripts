#!/bin/bash
set echo


#limit cputime 259200


echo `date`
matlab -nodisplay -nosplash <basinAve_fastpostProc.m  >ave_log
echo `date`


