#include "solution_io.h"
#include <algorithm>
#include <iostream>
#include <cassert>
int main(int argc, char** argv)
{
    if(argc == 1 or argc < 4 )
    {
        std::cerr << "Extract var fields from a SU2 dat file and save it to a solb/sol file" <<std::endl;
        std::cerr << "if varname is given extract specific variable otherwise the whole " <<
        "data minus the point coordinates is used " <<std::endl;
        std::cerr << "Usage " << argv[0] << " solution.dat mesh.su2 {sol,solb} [varname]"  <<std::endl;
        std::cerr << "Example " << argv[0] << " restart_solution.dat wing.su2 solb Mach"  <<std::endl;
        return 1 ;
    }


    std::vector<double> values; // flat array of all variables for all points
    std::vector<double> requested_variable_values;
    //std::map<std::string,int> varnames; // map from variables present in the file to column id
    std::vector<std::string> varnames; // map from variables present in the file to column id
    const char* filename = argv[1];
    const char* meshfile = argv[2];
    std::string varname;
    bool print_binary = true;

    std::string sarg(argv[3]);
    if( sarg == "sol")
    {
        print_binary = false;
    }
    else if(sarg != "solb")
    {
        std::cerr<< " Wrong output solution format specified, only sol,solb are supported \n";
        return 1;
    }


    if(argc == 5)
        varname = argv[4];
            
    /* strip ending and dot */
    std::string basefilename(filename); 
    basefilename.erase(basefilename.size() - 4, 4);

    read_variables_from_binary_restart_file(filename,meshfile,values,varnames);

    std::cout << "The solution file contains the following variables\n";
    for( auto item : varnames)
        std::cout << item << std::endl;

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
            std::size_t nPoints = values.size() / nvars;
            //std::cout << "nPoints " << nPoints << std::endl;
            requested_variable_values.resize(nPoints);
            basefilename.append("_");
            basefilename.append(varname);
            for(std::size_t i = 0 ; i < nPoints ; i++)
            {
                requested_variable_values[i] = values[ column + i*nvars];
                assert( requested_variable_values[i]  < 0.4);
            }

            if (print_binary) write_solb_with_scalar_vars(basefilename,1,requested_variable_values);
            else write_sol_with_scalar_vars(basefilename,1,requested_variable_values);
        }

    }
    else
    {
        int nvars = varnames.size() - 3;
        if (print_binary) write_solb_with_scalar_vars(basefilename,nvars,values);
        else write_sol_with_scalar_vars(basefilename,nvars,values);
    }

    return 0;
}