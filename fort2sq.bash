#!/bin/bash -f

nCorners=`head -1 fort.60`
echo $nCorners

inFile="dum.60"
outFile="sqgrid.in"

tail -n +2 fort.60 > $inFile

nx=`cat ./Include/gridparam.h | grep paramete |  cut -d "=" -f2 | cut -c1-3 `
ny=`cat ./Include/gridparam.h | grep paramete |  cut -d "=" -f3 | cut -c1-3 `

((nx++))
((ny++))

echo $nx $ny

if (( $nCorners > 7 )); then

	echo "doing the lat/lon edges thing"

	echo $ny >dum1.txt
	echo $nx >dum2.txt
	echo $ny >dum3.txt
	echo $nx >dum4.txt

	ii=1	
	while read line; do
		if (( $ii <= $ny )); then
			echo $line >> dum1.txt
#			echo $ii 'spacer' $line
#			echo $ii $ny
		fi

		if (( $ii > $ny - 1 && $ii <= $ny + $nx - 1 )); then
			echo $line >> dum2.txt
		fi

		if (( $ii > $ny + $nx - 2 && $ii <= $ny + $ny + $nx - 2 )); then
                	echo $line >> dum3.txt
        	fi

        	if (( $ii > $ny + $ny + $nx - 3 )); then
                	echo $line >> dum4.txt
       	 	fi

		((ii++))

	done <$inFile

	cat dum1.txt dum2.txt dum3.txt dum4.txt > sqgrid.in

	\rm dum1.txt dum2.txt dum3.txt dum4.txt $inFile

else

	echo "4 corners thing"

	# jgp this section works great for a simple 4-corner configuration 

	ii=1

	while read line; do
		if (( $ii == 1 )); then
			firstLine=`echo $line`
			echo "2"    > $outFile
			echo $line >> $outFile
		fi

		if (( $ii > 1 && $ii < 5 )); then
                	echo $line >> $outFile
			echo "2"   >> $outFile
                	echo $line >> $outFile
        	fi

        	((ii++))

	done <$inFile
        echo $firstLine >> $outFile

	\rm $inFile

fi
