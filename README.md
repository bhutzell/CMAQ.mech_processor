Test
# CMAQ.mech_processor
FORTRAN code for a utility, CHEMMECH, that creates RXNS Fortran90 modules for the CMAQ model version 5.1 as well as older code of the utility that creates RXDT.EXT and RXCM.EXT as in CMAQ version 4.7.1. 

This repository  contains the runscript called runit_chemmech and code directory to generates the RXNS_DATA_MODULE.F90 and RXNS_FUNC_MODULE.F90 files
for CMAQ version 5.1. THe repository also has a script called run.chemmech.CMAQv4.7.1 as well as inputs for using the older version of the CHEMMECH.

To use either version of the utility:
1) Compile the utility by modifying the Makefile in the BLD-kpp_hetero or BLD_chemmech_CMAQv4.7.1 directory. The changes set fortran compiler, cc compiler, their compilation
and link flags based on your system.

2) Next, modify the respective runscript for your own application input data or create test inputs the mechanism definitions files in sub directories under MECHS in the CMAQ.git repository or this repository. You may want to make a back up copy of the test data.
Also, define the output direct in the runscript. 

3) Execute the runscript and inspect the results under the output directory.
