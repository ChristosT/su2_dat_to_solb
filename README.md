# Convert between Binary SU2 restart files and INRIA/GMF .sol/solb file

## Build

    make all

## Usage

Extract all variables :

    ./dat2sol <filename.dat> <filename.sol/.solb> 

A file `metadata.txt` will be also written containing the names of each variable. If the sol/solb file is going to be converted again back todat the file is needed.

Extract a single variable :

    ./dat2sol <filename.dat> <filename.sol/.solb> <Variable name>

Convert solb back to SU2 restart file

    ./solb2dat <filename.solb> <mesh.su2> <newfile.dat> [ASCII]

* `mesh.su2` should match the `<filename.dat>` used to create the `solb` file. 
* `metadata.txt` from the previous step should be present.
Optionally passing `ASCII` will create an ascii solution file compatible with SU2 6.2.0 instead of a binary file.

## Notes
* Tested with SU2 6.2.0
