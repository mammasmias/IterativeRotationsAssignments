#! /bin/bash

# This script randomizes a given structure, and then runs the shape matching
# algorithms Arbalign, Fastoverlap, and IRA. 
# The argument passed to this script is the path to the xyz file containing 
# the original structure to be randomized and matched.
#
# Example run command:
#
#   ./run_single.sh data/lj_clusters/xyz/47.xyz 
#
# Before you run, make sure all software is compiled! SEE README

# the first argument passed is the path to structure
n=$1

#-----------------------------------
## Location of software files:
# Location of the ArbAlign software
aa_loc=./ArbAlign/

# Location of the fastoverlap software
fo_loc=./fastoverlap/fastoverlap-master/fastoverlap

# Location of the IRA software
ira_loc=../IRA/
#-----------------------------------



# copy the input file into temp_orig.xyz
cp ${n} temp_orig.xyz

# randomize the temp_orig.xyz, write into temp_random.xyz
./randomize/randomize.x < temp_orig.xyz > temp_random.xyz


#-----------------------------------
# Run the ArbAlign software:
#
# put both structures into principal axes for ArbAlign, this writes
# files temp_orig-prin.xyz and temp_random-prin.xyz, which contain the
# structures in their principal axes.
${aa_loc}/PrinCoords.py temp_orig.xyz
${aa_loc}/PrinCoords.py temp_random.xyz

#
# run ArbAlign with the structures in principal axes, write output to dump_aa
# This will also produce files reference.xyz and molecule-aligned-to-reference.xyz
${aa_loc}/ArbAlign-driver.py temp_orig-prin.xyz temp_random-prin.xyz > dump_aa


#-----------------------------------
# Run fastovelap software:
#
# convert from xyz to FO format, write temp_in.dat and temp_random.dat
./convert_xyz2FO/xyz2FO.x < temp_orig.xyz > temp_in.dat
./convert_xyz2FO/xyz2FO.x < temp_random.xyz > temp_random.dat

# run fastoverlap, this reads files temp_in.dat and temp_random.dat. 
# Write the output in dump_fo
python ${fo_loc}/sphericalAlignment_modread.py > dump_fo


#-----------------------------------
# Run IRA:
#
# Put the original and randomized xyz structures into the same file called temp_ira
cat temp_orig.xyz temp_random.xyz > temp_ira

# run IRA, write output dump_ira
${ira_loc}/ira_eq.x < temp_ira > dump_ira



#-----------------------------------
# Collect data from outputs:
rmsd_aa=$( grep "Best             RMSD" dump_aa | cut -c 23- )
rmsd_ira=$( grep rmsd dump_ira | cut -c 12- )
rmsd_fo=$( grep rmsd dump_fo | cut -c 11- )

#-----------------------------------
# write input file path, and values to screen
echo ${n} ${rmsd_aa} ${rmsd_fo} ${rmsd_ira}
#echo ${rmsd_ira}



