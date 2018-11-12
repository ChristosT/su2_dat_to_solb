#include "solution_io.h"
#include <iostream>
int main(int argc, char** argv)
{
    if(argc == 1 or argc < 4 )
    {
        std::cerr << "Convert a solb file holding scalar values of a solution prepared for SU2 to "
        << " native SU2 binary dat format";
        std::cerr << "Usage " << argv[0] << " solution.solb mesh.su2 newfile.dat [BINARY|ASCII]"  <<std::endl;
        return 1 ;
    }

    bool print_binary = true;
    if (argc == 5)
    {
        std::string sarg(argv[4]);
        if ( sarg == "ASCII")
            print_binary = false;
    }


    std::vector<double> values; // flat array of all variables for all points
    const char* solfile = argv[1];
    const char* meshfile = argv[2];
    const char* newfile = argv[3];

    int nvars;

    read_solb_with_scalar_vars(solfile,nvars,values);
    if(nvars + 3 != 18)
    {
        std::cerr << " Solution file should have 15 variables";
        return 1;
    }

    if(print_binary) write_variables_to_binary_restart_file(newfile,meshfile,values);
    else write_variables_to_ascii_restart_file(newfile,meshfile,values);

    return 0;
}
