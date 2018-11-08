#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#ifdef  NDEBUG

#define DPRINT(x) (void)(x);
#define SPRINT(x) (void)(x);

#else

#define DPRINT(x) std::cout << #x << " --> " << x << std::endl;
#define SPRINT(x) std::cout << #x << " --> |" << x <<"|" <<std::endl;

#endif
#define CHECK(x) if( not (x)) { std::cerr<< #x << " failed @" << __LINE__ << std::endl; std::abort();}
// read the number of points from an su2 mesh
std::size_t get_number_of_points(const char* filename)
{
    std::fstream file;
    file.open(filename,std::ios_base::in);
    std::string line;
    std::size_t num = 0;
    while(std::getline(file,line))
    {
        if ( line.size() > 0 )
        {
             int res = std::strncmp(line.c_str(),"NPOIN", std::min(line.size(),(long unsigned int)5));
             if (res == 0)
            {
                // split in spaces
                std::istringstream stream(line);
                stream >> line; // NPOINT=
                stream >> num;
                break;
            }
        }
    }
    file.close();

    return num;
}
// check CSolver::Read_SU2_Restart_Binary

/**
 * @brief read all variable value on every point from a su2 binary .dat file
 *
 * @param filename name of the file
 * @param nPoints point sin the mesh
 * @param[out] values  flat table of variable values in the form [ var0 for p0, var1, for p0,...]
 * @param[out] varnames a map of name of variables appearing in the file and their column ids
 */
void read_variables_from_binary_restart_file( const char* filename, std::size_t nPoints,
                                              std::vector<double>& values,
                                              std::map<std::string,int>& varnames)
{
    const int CGNS_STRING_SIZE = 33;
    std::FILE* pFile;
    pFile = std::fopen(filename,"rb");

    const int nRestart_Vars = 5;
    std::vector<int> Restart_Vars(nRestart_Vars); 

    // REad number of variables and make sure that it is 5
    std::size_t ret = std::fread(Restart_Vars.data(), sizeof(int), nRestart_Vars, pFile);

    CHECK( ret == nRestart_Vars);
    DPRINT(ret);

    // magic number in should be 535532
    CHECK(Restart_Vars[0] == 535532);
    DPRINT(Restart_Vars[0]);

    int  nFields = Restart_Vars[1];
    DPRINT(nFields);

    // Find the column where mach number is located
    char str_buf[CGNS_STRING_SIZE];
    for (int iVar = 0; iVar < nFields; iVar++) 
    {
        ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, pFile);
        CHECK(ret == CGNS_STRING_SIZE);
        varnames[str_buf] = iVar;
    }
    
    /*--- prepare return vector */

    values.resize(nFields*nPoints);

    /*--- Read in the data for the restart at all points. ---*/

    ret = std::fread(values.data(), sizeof(double), nFields*nPoints, pFile);
    CHECK(ret == nFields * nPoints);

    #ifndef NDEBUG
    //for( std::size_t i = 0 ; i < nPoints ; i++)
    for( std::size_t i = 0 ; i < 10 ; i++)
    {

        for( int j = 0 ; j < nFields ; j++)
        {
            std::cout << values[nFields*i + j] << " ";
        }
        std::cout << std::endl;
    }

    #endif
    
    /*--- Close the file. ---*/

    std::fclose(pFile);

}
void write_sol_with_scalar_vars(std::string filename,int nvars, std::vector<double>& values)
{
    std::fstream file;
    filename.append(".sol",4);
    file.open(filename,std::ios_base::out);
    file << "MeshVersionFormatted 2\n";
    file << "\n";
    file << "Dimension 3\n";
    file << "\n";
    file << "SolAtVertices \n" << values.size()/nvars << "\n";
    file << "1 ";
    for( int i = 0 ; i  < nvars; i++)
        file << "1 ";
    file << "\n";
    for(std::size_t i = 0 ; i < values.size(); )
    {
        for(int j = 0 ; j < nvars; j++ )
        {
            file << values[i] << "\n";
            i++;
        }
    }
    file << "\nEnd";
    file.close();
}

void write_solb_with_scalar_vars(std::string filename,int nvars, std::vector<double>& values)
{
    using std::FILE;
    using std::fwrite;
    using std::fopen;
    using std::fclose;
    
    FILE* pFile;
    filename.append(".solb",5);
    pFile = fopen(filename.c_str(),"wb");
    const int code = 1;
    const int version = 2;
    const int dimension_code = 3;
    const int dimension = 3;
    DPRINT(nvars);
    CHECK( not values.empty())
    int64_t res;

    // Write file header 
    {
        res = fwrite(&code, sizeof(int), 1, pFile);
        CHECK(res == 1);

        res = fwrite(&version, sizeof(int), 1, pFile);
        CHECK(res == 1);
        
        
        // end of the upcomming section
        int end_position = 3*sizeof(int) + ftell(pFile);
        fwrite(&dimension_code, sizeof(int), 1, pFile);

        fwrite(&end_position, sizeof(int), 1, pFile);

        fwrite(&dimension, sizeof(int), 1, pFile);
    }


    // Write solution header 

        const int keyword = 62 ; // = SolAtVertices
        const int npoints  = (int) values.size()/nvars ; // NOT SAFE for big meshes

        // end of the upcomming section
        int end_position = ftell(pFile) + (4 + nvars)*sizeof(int) + npoints*nvars*sizeof(double);

        fwrite(&keyword, sizeof(int), 1, pFile);
        fwrite(&end_position, sizeof(int), 1, pFile);


        fwrite(&npoints, sizeof(int), 1, pFile);

        const int solutions_per_node = nvars;
        fwrite(&solutions_per_node, sizeof(int), 1, pFile);

        const int solutions_type = 1; // 1 number per node
        for ( int i = 0 ; i < nvars ; i++)
            fwrite(&solutions_type, sizeof(int), 1, pFile);

    // write values at once
    res = fwrite(values.data(), sizeof(double)*values.size(), 1, pFile);
    CHECK(res == 1);

    // finalize file
    {
        const int keyword = 54 ; // =  GmfEmd
        fwrite(&keyword,sizeof(int),1,pFile);
        const int end_position = 0;
        fwrite(&end_position, sizeof(int), 1, pFile);
    }


    fclose(pFile);
}


int main(int argc, char** argv)
{
    if(argc == 1)
    {
        std::cerr << "Extract var fields from a SU2 dat file and save it ot a solb file" <<std::endl;
        std::cerr << "if varname is given extract specific variable" <<std::endl;
        std::cerr << "Usage " << argv[0] << "solution.dat mesh.su2 [varname]"  <<std::endl;
    }
    std::size_t nPoints = get_number_of_points(argv[2]);

    std::vector<double> values; // flat array of all variables for all points
    std::vector<double> requested_variable_values;
    std::map<std::string,int> varnames; // map from variables present in the file to column id
    const char* filename = argv[1];
    std::string varname;

    if(argc == 4)
        varname = argv[3];
            
    /* strip ending and dot */
    std::string basefilename(filename); 
    basefilename.erase(basefilename.size() - 4, 4);

    read_variables_from_binary_restart_file(filename,nPoints,values,varnames);
    DPRINT(values.size());

    // extract a single variable and print it in a sol/solb file
    if( not varname.empty())
    {
        auto iter = varnames.find(varname);
        if( iter == varnames.end())
        {
            std::cerr<< "This variables name does not exist in the su2 file" << std::endl;
            std::cerr<< "Variable names in provided su2 file are " << std::endl;
            for( auto item : varnames)
                std::cout << item.first << std::endl;
        }
        else
        {
            int column = iter->second;
            int nvars = varnames.size();
            requested_variable_values.resize(nPoints);
            for( int64_t i = 0 ; i < nPoints ; i++)
            {
                requested_variable_values[i] = values[ column + i*nvars];
            }

            write_sol_with_scalar_vars(basefilename,1,requested_variable_values);
            write_solb_with_scalar_vars(basefilename,1,requested_variable_values);
        }

    }
    else
    {
        int nvars = varnames.size();
        write_solb_with_scalar_vars(basefilename,nvars,values);
    }

    return 0;
}
