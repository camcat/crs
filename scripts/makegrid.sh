#/bin/bash

#Set ranges and resolution:
lat0=-43;
lat1=-41;
dlat=0.02;
lon0=172;
lon1=175.5;
dlon=0.02;
dep0=5;
dep1=25;
ddep=10;
Mmin=3;
Mmax=10;

for dep in $(seq $dep0 $ddep $dep1); do 
  for lon in $(seq $lon0 $dlon $lon1); do
    for lat in $(seq $lat0 $dlat $lat1); do
	 echo $lon $lat $dep $Mmin $Mmax
     done
  done
done

