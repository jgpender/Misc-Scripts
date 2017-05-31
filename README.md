# Misc-Scripts
Various scripts used for ROMS and MITgcm models

basinAve			basinAve_fastpostProc.m calculates basin average of 3D quantities
					temp, salt
				and surface quantities
					sst, sss, zeta, shflux, ssflux,swrad, sustr, svstr
				for each snapshot in a batch of ROMS output files.
				Plots of each quantity in ./figures.
	
				This can be done as ROMS is running. Start the script
					lurk.bash
				and it sleeps for an hour or two, wakes, runs basinAve_fastpostProc.m
				from the command line using bat.bash, start a new lurk job then quits.

basinAverages			defunct

checkModes*			This is work related to calculating the vertical mode structure from a 
				set of ROMs output.

Gridpack_chinook_template	Copy this when you want to make a new grid file	

Gridpack_template		defunct. This was for pacman.

gridTools			Simple script meant to help interalize ROMS terrain-following vertical grid.

HYCOMtemplate			Shell scripts for downloading HYCOM data.

matlab				An early version of the ROMS postprocessing tools. Up-to-date version is on Git.

postProcTools			Yet another version of the ROMS postprocessing tools on Git.

regridTools			Some scripts I haven't used in a while.

spectralTools			WIP scripts for analyzing vertical structure.

strat				Generates custom vertical structure for MITgcm input. gendata.m variant.

subregion			Script for subsetting some ROMS data

TPXOstuff			Scripts for doing tidal analysis of ROMS output. These are the main scripts
					tideComparison_K1.m
                                        tideComparison_M2.m
                                        tideComparison_S2.m
				and there are shell scripts you can use to run this from the command line.

