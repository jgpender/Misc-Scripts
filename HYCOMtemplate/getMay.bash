#!/bin/bash

if [ ! -d "data" ]; then 
mkdir data
fi

latMin=`grep latMin ./myLatLon.txt | cut -d "=" -f2` 
echo "latMin = " $latMin

lonMin=`grep lonMin ./myLatLon.txt | cut -d "=" -f2`
echo "lonMin = " $lonMin

latMax=`grep latMax ./myLatLon.txt | cut -d "=" -f2`
echo "latMax = " $latMax

lonMax=`grep lonMax ./myLatLon.txt | cut -d "=" -f2`
echo "lonMax = " $lonMax


part1='http://ncss.hycom.org/thredds/ncss/GLBa0.08/expt_91.2/'
# insert year
part2='?var=ssh&var=salinity&var=temperature&var=u&var=v'
part3='&north='$latMax'&west='$lonMin'&east='$lonMax'&south='$latMin
part4='&disableProjSubset=on&horizStride=1&time_start='
# insert date
part5='T00%3A00%3A00Z&time_end='
# insert date
part6='T00%3A00%3A00Z&timeStride=1&vertStride=1&addLatLon=true&accept=netcdf'


echo $part1
echo $part2
echo $part3
echo $part4
echo $part5

year=2016
monthStart='05'
dayStart=1
monthEnd='05'
dayEnd=28

#dayNumberStart=`date -d $year-$monthStart-$dayStart +%j`
#dayNumberEnd=`date -d $year-$monthEnd-$dayEnd +%j`
#echo $dayNumberStart
#echo $dayNumberEnd

echo $year
echo $monthStart
echo $monthEnd
echo $dayStart
echo $dayEnd

for i in `seq $dayStart $dayEnd`;
do
	if [  $i -lt 10 ] 
	then
		date="$year-$monthStart-0$i"
	else
		date="$year-$monthStart-$i"
	fi

	echo $date

        nDay=`date -d $date +%j`

        outFile="./data/HYCOM_GLBa0.08_"$year"_"$nDay".nc"
	echo $outFile

	myURL=$part1$year$part2$part3$part4$date$part5$date$part6

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

done
