#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
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
void read_variable_from_binary_restart_file( const char* filename, const char* fieldname, std::vector<double>& mach)
{
    const int CGNS_STRING_SIZE = 33;
    std::FILE* pFile;
    pFile = std::fopen(filename,"rb");
    std::size_t nPoints = mach.size();
    int len = std::min((int)std::strlen(fieldname),CGNS_STRING_SIZE);
    DPRINT(len);
    DPRINT(fieldname);

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
    int column = -1;
    for (int iVar = 0; iVar < nFields; iVar++) 
    {
        ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, pFile);
        CHECK(ret == CGNS_STRING_SIZE);
        if ( std::strncmp(fieldname,str_buf,len) == 0)
            column = iVar;

        //DPRINT(str_buf);
        DPRINT("-------------");
    }
    CHECK(column != -1);
    DPRINT(column);
    /*--- For now, create a temp 1D buffer to read the data from file. ---*/

    std::vector<double> Restart_Data(nFields*nPoints);

    /*--- Read in the data for the restart at all local points. ---*/

    ret = std::fread(Restart_Data.data(), sizeof(double), nFields*nPoints, pFile);
    CHECK(ret == nFields * nPoints);

    #ifndef NDEBUG
    //for( std::size_t i = 0 ; i < nPoints ; i++)
    for( std::size_t i = 0 ; i < 10 ; i++)
    {

        for( int j = 0 ; j < nFields ; j++)
        {
            std::cout << Restart_Data[nFields*i + j] << " ";
        }
        std::cout << std::endl;
    }

    #endif
    
    for( std::size_t i = 0 ; i < nPoints ; i++)
        mach[i] = Restart_Data[nFields*i + column];

    /*--- Close the file. ---*/

    std::fclose(pFile);

}
void write_sol_with_var(std::string filename,std::vector<double>& values)
{
    std::fstream file;
    filename.append(".sol",4);
    file.open(filename,std::ios_base::out);
    file << "MeshVersionFormatted 2\n";
    file << "\n";
    file << "Dimension 3\n";
    file << "\n";
    file << "SolAtVertices \n" << values.size() << "\n";
    file << "1 1\n";
    for(double x : values)
        file << x << "\n";
    file << "\nEnd";
    file.close();
}
void write_solb_with_var(std::string filename,std::vector<double>& values)
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

    // Write file header 
    {
        fwrite(&code, sizeof(int), 1, pFile);
        fwrite(&version, sizeof(int), 1, pFile);
        
        
        // end of the upcomming section
        int end_position = 3*sizeof(int) + ftell(pFile);
        fwrite(&dimension_code, sizeof(int), 1, pFile);

        fwrite(&end_position, sizeof(int), 1, pFile);

        fwrite(&dimension, sizeof(int), 1, pFile);
    }


    // Write solution header 

        const int keyword = 62 ; // = SolAtVertices
        const int npoints  = (int) values.size(); // NOT SAFE for big meshes

        // end of the upcomming section
        int end_position = ftell(pFile) + 5*sizeof(int) + npoints* sizeof(double);

        fwrite(&keyword, sizeof(int), 1, pFile);
        fwrite(&end_position, sizeof(int), 1, pFile);
        fwrite(&npoints, sizeof(int), 1, pFile);

        const int solutions_per_node = 1;
        fwrite(&solutions_per_node, sizeof(int), 1, pFile);

        const int solutions_type = 1; // 1 number per node
        fwrite(&solutions_type, sizeof(int), 1, pFile);

    // write values at once
    fwrite(values.data(), sizeof(double)*npoints, 1, pFile);

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
        std::cerr << "Extract mach field from a SU2 dat file and save it ot a solb file" <<std::endl;
        std::cerr << "Usage " << argv[0] << "solution.dat mesh.su2"  <<std::endl;
    }
    std::size_t nPoints = get_number_of_points(argv[2]);

    std::vector<double> mach(nPoints);
    const char* filename = argv[1];

    read_variable_from_binary_restart_file(filename,"Mach",mach);

    std::string basefilename(filename);
    //strip .dat
    basefilename.erase(basefilename.size() - 4, 4);

    write_sol_with_var(basefilename,mach);
    write_solb_with_var(basefilename,mach);

    return 0;
}
