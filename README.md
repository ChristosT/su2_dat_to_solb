# dat2solb

Convert between Binary SU2 restart files and INRIA/GMF .sol/solb file

## Build

    make all

## Usage

Extract all variables :

    ./dat2sol <filename.dat> <filename.sol/.solb> 

Two temporary files will be also written:

* `varnames.txt` containing the name of each variable. If the sol/solb file is going to be converted again back this file is needed.
* `metadata.txt` in case the .dat file contains aditional information like angle of attack, exit iteration etc (default in SU2 6.2.0) this file will hold its values. If the sol/solb file is going to be converted again back this file is needed.

Extract a single variable :

    ./dat2sol <filename.dat> <filename.sol/.solb> <Variable name>

Convert solb (holding all variables) back to SU2 restart file

    ./solb2dat <filename.solb> <mesh.su2> <newfile.dat> [ASCII]

* `mesh.su2` should match the `<filename.dat>` used to create the `solb` file. 
* If `metadata.txt` was generated in previous step it should be present.
Optionally passing `ASCII` will create an ascii solution file compatible with SU2 6.2.0 instead of a binary file.

## Notes
* Tested with SU2 6.2.0, 7.0.1
