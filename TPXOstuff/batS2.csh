#!/bin/csh
source ~/.runROMS
set echo

limit cputime 259200

pwd


matlab -nodisplay -nosplash <  tideComparison_S2.m
