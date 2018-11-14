#include "solution_io.h"
#include <algorithm>
#include <iostream>
#include <cassert>
int main(int argc, char** argv)
{
    if(argc == 1 or argc < 3 )
    {
        std::cerr << "Extract var fields from a SU2 dat file and save it to a solb/sol file" <<std::endl;
        std::cerr << "if varname is given extract specific variable otherwise the whole " <<
        "data minus the point coordinates is used " <<std::endl;
        std::cerr << "Usage " << argv[0] << " solution.dat outfile.{sol,solb} [varname]"  <<std::endl;
        std::cerr << "Example " << argv[0] << " restart_solution.dat restart_solution.solb Mach"  <<std::endl;
        return 1 ;
    }


    std::vector<double> values; // flat array of all variables for all points
    std::vector<double> requested_variable_values;
    //std::map<std::string,int> varnames; // map from variables present in the file to column id
    std::vector<std::string> varnames; // map from variables present in the file to column id
    const char* filename = argv[1];
    std::string varname;
    std::size_t nPoints;
    bool print_binary = true;

    std::string outfile(argv[2]);
    int pos = outfile.find_last_of(".");
    std::string suffix = outfile.substr(pos);
    if( suffix == ".sol")
    {
        print_binary = false;
    }
    else if(suffix != ".solb")
    {
        std::cerr<< " Wrong output solution format specified, only sol,solb are supported \n";
        return 1;
    }


    if(argc == 4)
        varname = argv[3];
            
    /* strip ending and dot */
    std::string basefilename(filename); 
    basefilename.erase(basefilename.size() - 4, 4);

    read_variables_from_binary_restart_file(filename,values,varnames,nPoints);

    std::cout << "The solution file contains the following variables\n";
    std::cout << "--------------------------\n";
    for( auto item : varnames)
        std::cout << item << std::endl;
    std::cout << "--------------------------\n";

    // extract a single variable and print it in a sol/solb file
    if( not varname.empty())
    {
        auto iter = std::find(varnames.begin(),varnames.end(),varname);
        if( iter == varnames.end())
        {
            std::cerr<< "This variables name does not exist in the su2 file" << std::endl;
            std::cerr<< "Variable names in provided su2 file are " << std::endl;
            for( auto item : varnames)
                std::cout << item << std::endl;
        }
        else
        {
            int column = iter - varnames.begin() - 3; // -3 because coornates are disgraded
            if( column  < 0 )
            {
                std::cerr<< "coordinates are disregarded during reading " << std::endl;
                return 1;
            }
            std::cout << "column " << column << std::endl;
            int nvars = varnames.size() - 3;
            std::cout << "nPoints " << nPoints << std::endl;
            std::cout << "nvars " << nvars << std::endl;
            requested_variable_values.resize(nPoints);

            for(std::size_t i = 0 ; i < nPoints ; i++)
            {
                requested_variable_values[i] = values[ column + i*nvars];
            }

            if (print_binary) write_solb_with_scalar_vars(outfile,1,nPoints,requested_variable_values);
            else write_sol_with_scalar_vars(outfile,1,nPoints, requested_variable_values);
        }

    }
    else
    {
        int nvars = varnames.size() - 3;
        if (print_binary) write_solb_with_scalar_vars(outfile,nvars,nPoints,values);
        else write_sol_with_scalar_vars(outfile,nvars,nPoints,values);
    }

    return 0;
}
