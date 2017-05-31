#!/bin/bash

if [ ! -d "data" ]; then 
mkdir data
fi

latMin=`grep latMin ./myLatLon.txt | cut -d "=" -f2` 
latMin=`echo "$latMin - .08"  |bc`
echo "latMin = " $latMin

lonMin=`grep lonMin ./myLatLon.txt | cut -d "=" -f2`
lonMin=`echo "$lonMin - .08"  |bc`
echo "lonMin = " $lonMin

latMax=`grep latMax ./myLatLon.txt | cut -d "=" -f2`
latMax=`echo "$latMax + .08"  |bc`
echo "latMax = " $latMax

lonMax=`grep lonMax ./myLatLon.txt | cut -d "=" -f2`
lonMax=`echo "$lonMax + .08"  |bc`
echo "lonMax = " $lonMax


#part1='http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.2?var=surf_el&var=water_temp&north='$latMax'&west='$lonMin'&east='$lonMax'&south='$latMin'&disableProjSubset=on&horizStride=1&time_start='
#part2='T00%3A00%3A00Z&time_end='
#part3='T00%3A00%3A00Z&timeStride=1&vertCoord=&accept=netcdf'

part1='http://ncss.hycom.org/thredds/ncss/GLBa0.08/expt_91.2/2016?var=ssh&var=salinity&var=temperature&var=u&var=v'
part2='&north='$latMax'&west='$lonMin'&east='$lonMax'&south='$latMin
part3='&disableProjSubset=on&horizStride=1&time_start='
part4='T00%3A00%3A00Z&time_end='
part5='T00%3A00%3A00Z&timeStride=1&vertStride=1&addLatLon=true&accept=netcdf'

echo $part1
echo $part2
echo $part3
echo $part4
echo $part5


# pick any file for the "grid" file
year=2016
month='04'
day='25'

echo $year
echo $month
echo $day

date="$year-$month-$day"
echo $date


outFile="HYCOM_GLBa0.08_PALAU_grid.nc"
echo $outFile

myURL=$part1$part2$part3$date$part4$date$part5

echo $myURL


        wget -O $outFile $myURL

        ncrename -O -h -d MT,ocean_time         $outFile
        ncrename -O -h -d Depth,z               $outFile

        ncrename -O -h -v MT,ocean_time         $outFile
        ncrename -O -h -v Depth,z               $outFile

     	ncrename -O -h -v temperature,temp   	$outFile
        ncrename -O -h -v salinity,salt		$outFile

        ncrename -O -h -v Latitude,lat        	$outFile
        ncrename -O -h -v Longitude,lon       	$outFile
