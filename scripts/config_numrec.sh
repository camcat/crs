#/bin/bash

# This script sets the folder containing source code for Numerical Recipes, or allows to compile the code without using them.
# You only need the Numerical Recipes if you want to use the option 'splines' for aseismic stress sources. If you are not including aseismic processes, or if you use another option ('log' or 'lin' in the parameter file), you can set: "numerical_recipes_installed=no" and run this script. Otherwise, you need to set the paths to the Numerical Recipes source files below.
# We hope to make the code completely Numerical-Recipes free in the future!
#
# The variable numrec_folder is the path to your numerical recipes source files (the root folder, without slash at the end).
# 
# Notes: 
# 1. you need to run this scripts from the root crs folder.
# 2. For the changes to take effect, after running the script clean and recompile the code, as follows:
#	cd Release; make clean; make all; cd ..
#

# Write yes or no here:
numerical_recipes_installed=no

# Put the path to the Numerical Recipes source code here:
numrec_folder="/home/camilla/code/C/NR"


#---------------Don't need to change anything below here--------------------------------%

mainfile="Release/code/subdir.mk"
utilfile="Release/code/src/util/subdir.mk"
numrecfile="Release/code/src/nr/subdir.mk"
# Name of flag to deactivae numerical recipes.
noNRflag="_no_numericalrecipes"

if [ $numerical_recipes_installed = "yes" ];
  then
	#Remove flag which deactivates numerical recipes:
        grep -rl "$noNRflag" $utilfile | xargs sed -i "s* -D$noNRflag**g"
        grep -rl "$noNRflag" $mainfile | xargs sed -i "s* -D$noNRflag**g"
        echo 
	echo " Removed flag $noNRflag from $utilfile, $mainfile."

	#Comment out previous paths in numerical recipes makefile:
	sed -i "s*^numrec_path*#numrec_path*g" $numrecfile
	#Add a line with the path to numerical recipes to the local numerical recipes makefile:
	sed -i "5s*^*numrec_path=$numrec_folder\n*" $numrecfile

	echo " The following include pattern has been set in $numrecfile: $numrec_folder"
	echo
	
elif [ $numerical_recipes_installed = "no" ];
  then

	#determine which compiler is set (gcc or mpicc)

	find_gcc=$(grep -rl gcc $utilfile | grep 'subdir\|makefile' | wc -l)
	find_mpicc=$(grep -rl mpicc $utilfile | grep 'subdir\|makefile' | wc -l)

	if [ $find_gcc != "0" ] && [ $find_mpicc != "0" ]
	then
	   echo 
	   echo " Could not determine which compiler is being used: please set it manually in config_numrec.sh"
	   echo
	else
	   sed -i "s/^numrec_path/#numrec_path/g" $numrecfile
	   if [ $find_gcc != "0" ]
	   then
		sed -i "s*gcc*gcc -D$noNRflag*" $utilfile
		sed -i "s*gcc*gcc -D$noNRflag*" $mainfile
	   else
		sed -i "s*mpicc*mpicc -D$noNRflag*" $utilfile
		sed -i "s*mpicc*mpicc -D$noNRflag*" $mainfile
	   fi
	   echo
	   echo " Commented out numrec_path in $numrecfile, $mainfile"
	   echo
	   echo " Added flag: $noNRflag to $utilfile."
	   echo " Now you can compile without the Numerical Recipes source code."
	   echo 	
	fi	

else
  echo
  echo " Illegal value for variable 'numerical_recipes_installed' in config_numrec.sh. Should be 'yes' or 'no' "
  echo
fi   
  
