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
#include <cassert>
#ifdef  NDEBUG

#define DPRINT(x) (void)(x);
#define SPRINT(x) (void)(x);

#else

#define DPRINT(x) std::cout << #x << " --> " << x << std::endl;
#define SPRINT(x) std::cout << #x << " --> |" << x <<"|" <<std::endl;

#endif
#define CHECK(x) if( not (x)) { std::cerr<< #x << " failed @" << __LINE__ << std::endl; std::abort();}

struct MetaData
{
    /*--- External iteration ---*/ 
    // EXT_ITER
    int ext_iter;
    /*--- Angle of attack ---*/
    // AOA
    double aoa;
    /*--- Sideslip angle ---*/
    // SIDESLIP_ANGLE
    double sangle;
    /*--- BCThrust angle ---*/
    // INITIAL_BCTHRUST
    double initial_bc_thrust;

    double other[5];
    void pack()
    {
        std::ofstream file;
        file.open("metadata.txt");
        CHECK(file.is_open());
        file << ext_iter << " ";
        file << aoa << " ";
        file << sangle << " ";
        file << initial_bc_thrust << " ";
        for( int i = 0 ; i < 5 ; i++)
            file << other[i] << " ";
        file.close();
    }
    void unpack()
    {
        std::ifstream file;
        file.open("metadata.txt");
        CHECK(file.is_open());
        file >> ext_iter;
        file >> aoa;
        file >> sangle;
        file >> initial_bc_thrust;
        for( int i = 0 ; i < 5 ; i++)
        file >> other[i];
        file.close();
    }

};

// read the number of points from an su2 mesh
static std::size_t get_number_of_points(const char* filename)
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

static void pack_variable_names(std::vector<std::string>& varnames)
{
    std::ofstream file;
    file.open("varnames.txt");
    CHECK(file.is_open());
    file << varnames.size() << " ";
    for(std::string& name :varnames)
        file << name << " ";

    file.close();
}
static void unpack_variable_names(std::vector<std::string>& varnames)
{
    std::ifstream file;
    file.open("varnames.txt");
    CHECK(file.is_open());

    std::size_t size;
    file >> size;
    varnames.resize(size);
    for(std::size_t i = 0 ; i < size; i++)
        file >> varnames[i];

    file.close();
}


// check CSolver::Read_SU2_Restart_Binary in SU2 for more

void read_variables_from_binary_restart_file( const char* filename, 
                                              const char* meshfile,
                                              std::vector<double>& values,
                                              std::vector<std::string>& varnames,
                                              std::size_t& nPoints)
{
    nPoints = get_number_of_points(meshfile);
    DPRINT(nPoints);
    const int CGNS_STRING_SIZE = 33;
    std::FILE* pFile;
    pFile = std::fopen(filename,"rb");

    const int nRestart_Vars = 5;
    std::vector<int> Restart_Vars(nRestart_Vars); 

    // REad number of variables and make sure that it is 5
    std::size_t ret = std::fread(Restart_Vars.data(), sizeof(int), nRestart_Vars, pFile);

    CHECK( ret == nRestart_Vars);
    DPRINT(ret);

    // magic number of SU2 file should be 535532
    CHECK(Restart_Vars[0] == 535532);

    int  nFields = Restart_Vars[1];
    DPRINT(nFields);
    varnames.resize(nFields);
    int nVars = nFields - 3; // disregard x,y,z coordinates

    // Find the column where mach number is located
    char str_buf[CGNS_STRING_SIZE];
    for (int iVar = 0; iVar < nFields; iVar++) 
    {
        ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, pFile);
        CHECK(ret == CGNS_STRING_SIZE);
        varnames[iVar] = str_buf;
    }
    
    /*--- prepare return vector */

    values.resize(nVars*nPoints);

    /*--- Read in the data for the restart at all points. ---*/

    bool skip_coords = true; // if we need all data we can make it a flag
    if( not skip_coords)
    {
        ret = fread(values.data(), sizeof(double), nFields*nPoints, pFile);
        CHECK(ret == nFields * nPoints);
    }
    else
    {
        double coords[3];
        for( size_t i = 0 ; i < nPoints;i++)
        {
            ret = fread(coords, sizeof(double), 3, pFile);
            CHECK(ret == 3);
            ret = fread(&values[i*nVars], sizeof(double)*nVars, 1, pFile);
            CHECK(ret == 1);
        }
    }
// https://github.com/su2code/SU2/blob/5d5571f7fc9e4b0f77d093b31d593e5fc94f426f/SU2_CFD/src/solver_structure.cpp#L2608
    // read Metadata
    fseek(pFile, -(sizeof(int) + 8 * sizeof(double)),SEEK_END);
    MetaData md;
    
    ret = fread(&md.ext_iter, sizeof(int), 1, pFile);
    CHECK( ret == 1)
    ret = fread(&md.aoa, sizeof(double), 1, pFile);
    CHECK( ret == 1)
    ret = fread(&md.sangle, sizeof(double), 1, pFile);
    CHECK( ret == 1)
    ret = fread(&md.initial_bc_thrust, sizeof(double), 1, pFile);
    CHECK( ret == 1)
    ret = fread(&md.other, 5*sizeof(double), 1, pFile);
    CHECK( ret == 1)
    md.pack();

    #ifndef NDEBUG
    //for( std::size_t i = 0 ; i < nPoints ; i++)
    //for( std::size_t i = 0 ; i < 10 ; i++)
    //{

    //    for( int j = 0 ; j < nVars ; j++)
    //    {
    //        printf("%f ",values[nVars*i + j]);
    //    }
    //    std::cout << std::endl;
    //}

    #endif
    
    /*--- Close the file. ---*/

    std::fclose(pFile);

    pack_variable_names(varnames);

}

void write_variables_to_binary_restart_file( const char* solutionfile, 
                                             const char* meshfile,
                                             std::vector<double>& values)
{
    using std::FILE;
    using std::fopen;
    using std::fwrite;
    
    const int CGNS_STRING_SIZE = 33;

    // get points from mesh file
    std::fstream mesh;
    mesh.open(meshfile);

    std::string line;

    int ret;
    CHECK(mesh.is_open());
    DPRINT(meshfile);
    while(std::getline(mesh,line))
    {
        // Check whether the line starts by a character
        // of if it is a comment
        if(line.size() > 0 and (isalpha(line[0])) )
        {
            ret = std::strncmp(line.c_str(),"NPOIN", std::min((int)line.size(),5) );
            if( ret == 0 )
            {

                break;
            }
        }
    }
    CHECK(ret == 0);
    std::size_t nPoints = 0;
    std::vector< std::array<double,3> > points;
    std::stringstream ss(line);
    ss >> line; // NPOIN=
    ss >> nPoints; 
    DPRINT(nPoints);
    double x,y,z;
    points.reserve(nPoints);
    while(std::getline(mesh,line))
    {
        std::stringstream ss(line);
        ss >> x >> y >> z;

        points.push_back({x,y,z});
    }

    mesh.close();

    // write the solution file
    FILE* pFile;
    pFile = fopen(solutionfile,"wb");



    int nVars = values.size() / nPoints;
    const int nRestart_Vars = 5;
    int Restart_Vars[nRestart_Vars] = {-1};
    int nFields = 3 + nVars;
    Restart_Vars[0] = 535532 ; // SU2 magic numbers
    Restart_Vars[1] = nFields;
    
    ret = fwrite(Restart_Vars, sizeof(int), nRestart_Vars, pFile);
    CHECK( ret == nRestart_Vars);

    CHECK( nVars + 3 == 18);
    std::vector<std::string> varnames;
    unpack_variable_names(varnames);

    char fieldname[CGNS_STRING_SIZE];
        
    for (int iVar = 0; iVar < nFields; iVar++) 
    {
        std::strncpy(fieldname, varnames[iVar].data(), std::min((int)varnames.size(),CGNS_STRING_SIZE));
        ret = fwrite(fieldname, sizeof(char), CGNS_STRING_SIZE, pFile);
        CHECK(ret == CGNS_STRING_SIZE);
    }
    
    
    for( std::size_t i = 0 ; i < nPoints ; i++)
    {
        fwrite(points[i].data(), sizeof(double), 3, pFile);
        fwrite(&values[i*nVars], sizeof(double),nVars,pFile);
    }
    MetaData md;
    md.unpack();
    DPRINT(md.ext_iter);
    ret = fwrite(&md.ext_iter, sizeof(int), 1, pFile);
    CHECK( ret == 1)
    ret = fwrite(&md.aoa, sizeof(double), 1, pFile);
    CHECK( ret == 1)
    ret = fwrite(&md.sangle, sizeof(double), 1, pFile);
    CHECK( ret == 1)
    ret = fwrite(&md.initial_bc_thrust, sizeof(double), 1, pFile);
    CHECK( ret == 1)
    ret = fwrite(md.other, sizeof(double)*5, 1, pFile);
    CHECK( ret == 1)


    /*--- Close the file. ---*/

    fclose(pFile);

}

void write_sol_with_scalar_vars(std::string filename,int nVars, std::size_t nPoints, std::vector<double>& values)
{
    std::fstream file;
    filename.append(".sol",4);
    file.open(filename,std::ios_base::out);
    file << "MeshVersionFormatted 2\n";
    file << "\n";
    file << "Dimension 3\n";
    file << "\n";
    file << "SolAtVertices \n" << nPoints << "\n";
    file << "1 ";
    for( int i = 0 ; i  < nVars; i++)
        file << "1 ";
    file << "\n";
    for(std::size_t i = 0 ; i < values.size(); )
    {
        for(int j = 0 ; j < nVars; j++ )
        {
            file << values[i] ;
            i++;
        }
        file << "\n";
    }
    file << "\nEnd";
    file.close();
}

void write_solb_with_scalar_vars(std::string filename,int nVars, std::size_t nPoints, std::vector<double>& values)
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
    DPRINT(nVars);
    DPRINT(nPoints);
    DPRINT(values.size());
    CHECK( not values.empty())
    CHECK(values.size()  == nVars *nPoints);
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

        // end of the upcomming section
        int end_position = ftell(pFile) + (4 + nVars)*sizeof(int) + nPoints*nVars*sizeof(double);

        fwrite(&keyword, sizeof(int), 1, pFile);
        fwrite(&end_position, sizeof(int), 1, pFile);


        int nPoints_32 = (int) nPoints; //TODO  can we use big integers ??
        DPRINT(nPoints_32);
        fwrite(&nPoints_32, sizeof(int), 1, pFile);

        const int solutions_per_node = nVars;
        fwrite(&solutions_per_node, sizeof(int), 1, pFile);

        const int solutions_type = 1; // 1 number per node
        for ( int i = 0 ; i < nVars ; i++)
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
void read_solb_with_scalar_vars(std::string filename,int& nVars, std::vector<double>& values)
{
    using std::FILE;
    using std::fwrite;
    using std::fopen;
    using std::fclose;
    
    FILE* pFile;
    pFile = fopen(filename.c_str(),"rb");
    int code;
    int version;
    int dimension_code;
    int dimension;
    CHECK(values.empty())
    int64_t res;

    // Write file header 
    {
        res = fread(&code, sizeof(int), 1, pFile);
        CHECK(res == 1);
        CHECK(code == 1);

        res = fread(&version, sizeof(int), 1, pFile);
        CHECK(res == 1);
        CHECK(version == 2 || version == 3);
        
        
        // end of the upcomming section
        int end_position;
        fread(&dimension_code, sizeof(int), 1, pFile);
        CHECK(dimension_code == 3);

        fread(&end_position, sizeof(int), 1, pFile);
        CHECK(end_position > 0);

        fread(&dimension, sizeof(int), 1, pFile);
        CHECK(dimension == 3);

        // size of next
        //end_position -= 3*sizeof(int);
    }


    // Write solution header 

        int keyword;
        int nPoints;
        int end_position;


        fread(&keyword, sizeof(int), 1, pFile);
        CHECK(keyword == 62) ; // SolAtVertices = 62

        // end of the upcomming section
        fread(&end_position, sizeof(int), 1, pFile);


        fread(&nPoints, sizeof(int), 1, pFile);
        CHECK(nPoints > 0);

        int solutions_per_node;
        fread(&solutions_per_node, sizeof(int), 1, pFile);
        CHECK(solutions_per_node > 0);

        nVars = solutions_per_node;

        //int end_position = ftell(pFile) + (4 + nvars)*sizeof(int) + npoints*nvars*sizeof(double);
        int solutions_type; 
        for ( int i = 0 ; i < nVars ; i++)
        {
            fread(&solutions_type, sizeof(int), 1, pFile);
            CHECK( solutions_type == 1);
        }

        values.resize( nVars * nPoints);


        // read values at once
        res = fread(values.data(), sizeof(double)*values.size(), 1, pFile);
        CHECK(res == 1);

    // finalize file
    {
        int keyword;
        fread(&keyword,sizeof(int),1,pFile);
        CHECK( keyword == 54) ; // = GmfEmd

        int end_position;
        fread(&end_position, sizeof(int), 1, pFile);
        CHECK(end_position == 0);
    }


    fclose(pFile);
}
