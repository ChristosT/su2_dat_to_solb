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
#include <iomanip>
#ifdef  NDEBUG

#define DPRINT(x) (void)(x);
#define SPRINT(x) (void)(x);

#else

#define DPRINT(x) std::cout << #x << " --> " << x << std::endl;
#define SPRINT(x) std::cout << #x << " --> |" << x <<"|" <<std::endl;

#endif
#define CHECK(x) if( not (x)) { std::cerr<< #x << " failed @" << __FILE__ << ":" << __LINE__ << std::endl; std::abort();}

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
    bool unpack()
    {
        std::ifstream file;
        file.open("metadata.txt");
        if( not file.good())
            return false;
        CHECK(file.is_open());
        file >> ext_iter;
        file >> aoa;
        file >> sangle;
        file >> initial_bc_thrust;
        for( int i = 0 ; i < 5 ; i++)
            file >> other[i];
        file.close();
        return true;
    }

};
namespace 
{
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

    void pack_variable_names(std::vector<std::string>& varnames)
    {
        std::ofstream file;
        file.open("varnames.txt");
        CHECK(file.is_open());
        file << varnames.size() << " ";
        for(std::string& name :varnames)
            file << name << " ";

        file.close();
    }
    void unpack_variable_names(std::vector<std::string>& varnames)
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

    int determine_dimension(std::vector<std::string>& varnames)
    {
        for(std::string& name :varnames)
            if (name == "z" || name == "Z")
                return 3;
        return 2;
    }
}


// check CSolver::Read_SU2_Restart_Binary in SU2 for more

void read_variables_from_binary_restart_file( const char* filename, 
                                              std::vector<double>& values,
                                              std::vector<std::string>& varnames,
                                              std::size_t& nPoints,
                                              int& dimension)
{
    using std::FILE;
    using std::fopen;
    using std::fclose;
    using std::fwrite;

    const int CGNS_STRING_SIZE = 33;
    FILE* pFile;
    pFile = fopen(filename,"rb");
    if(pFile == NULL)
    {
        std::cerr << __func__ << "@" << __LINE__ << " cannot open file " << filename  << std::endl;
        std::abort();
    }

    const int nRestart_Vars = 5;
    std::vector<int> Restart_Vars(nRestart_Vars); 

    /* SU2 binary file header has the following variables in the header
     * [ SU2 magic number, number of variables per point, number of points, metatadata_ints, metatadata_doubles ]
     * - SU2 magic number is expected to be 535532
     * - metatadata_ints is expected to be  1
     * - metatadata_doubles is expected to be 8
     *
     * for more see
     * see void COutput::WriteRestart_Parallel_Binary
     *  at https://github.com/su2code/SU2/blob/master/SU2_CFD/src/output_structure.cpp#L17527
     */

    // Read number of variables and make sure that it is 5
    std::size_t ret = fread(Restart_Vars.data(), sizeof(int), nRestart_Vars, pFile);
    CHECK( ret == nRestart_Vars);
    DPRINT(ret);

    // magic number of SU2 file should be 535532
    CHECK(Restart_Vars[0] == 535532);

    int  nFields = Restart_Vars[1];
    nPoints = Restart_Vars[2];

    DPRINT(nFields);
    varnames.resize(nFields);

    // Find the column where mach number is located
    char str_buf[CGNS_STRING_SIZE];
    for (int iVar = 0; iVar < nFields; iVar++) 
    {
        ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, pFile);
        CHECK(ret == CGNS_STRING_SIZE);
        varnames[iVar] = str_buf;
    }
    dimension = determine_dimension(varnames);
    int nVars = nFields - dimension; // disregard x,y,z coordinates

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
            ret = fread(coords, sizeof(double), dimension, pFile);
            CHECK(ret == (std::size_t)dimension);
            ret = fread(&values[i*nVars], sizeof(double)*nVars, 1, pFile);
            CHECK(ret == 1);
        }
    }

    // read Metadata
    // see also https://github.com/su2code/SU2/blob/5d5571f7fc9e4b0f77d093b31d593e5fc94f426f/SU2_CFD/src/solver_structure.cpp#L2608
    int metatadata_ints = Restart_Vars[3];
    int metatadata_doubles = Restart_Vars[4];


    if(metatadata_ints > 0)
    {
        CHECK(metatadata_ints == 1);
        fseek(pFile, -(metatadata_ints * sizeof(int) + metatadata_doubles * sizeof(double)),SEEK_END);
        MetaData md;

        ret = fread(&md.ext_iter, sizeof(int), 1, pFile);
        CHECK( ret == 1);
        if(metatadata_doubles > 0 )
        {
            CHECK(metatadata_doubles == 8);
            ret = fread(&md.aoa, sizeof(double), 1, pFile);
            CHECK( ret == 1);
            ret = fread(&md.sangle, sizeof(double), 1, pFile);
            CHECK( ret == 1);
            ret = fread(&md.initial_bc_thrust, sizeof(double), 1, pFile);
            CHECK( ret == 1);
            ret = fread(&md.other, 5*sizeof(double), 1, pFile);
            CHECK( ret == 1);
        }
        md.pack();
    }

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
    int dimension = 0;
    std::size_t nPoints = 0;

    // get points from mesh file
    std::fstream mesh;
    mesh.open(meshfile);

    std::string line;

    int ret;
    CHECK(mesh.is_open());
    DPRINT(meshfile);
    while(std::getline(mesh,line) && (nPoints == 0 || dimension == 0))
    {
        // Check whether the line starts by a character
        // or if it is a comment
        if(line.size() > 0 and (isalpha(line[0])) )
        {
            ret = std::strncmp(line.c_str(),"NDIME", std::min((int)line.size(),5) );
            if( ret == 0 )
            {
                std::stringstream ss(line);
                ss >> line; // NDIME=
                ss >> dimension;
            }
            ret = std::strncmp(line.c_str(),"NPOIN", std::min((int)line.size(),5) );
            if( ret == 0 )
            {
                std::stringstream ss(line);
                ss >> line; // NPOIN=
                ss >> nPoints;
            }
        }
    }
    CHECK(ret == 0);
    std::vector< double > points;
    DPRINT(nPoints);
    double x;
    points.reserve(nPoints*dimension);
    while(std::getline(mesh,line))
    {
        std::stringstream ss(line);
        for(int i = 0 ; i < dimension; i++)
        {
            ss >> x;
            points.push_back(x);
        }
    }

    mesh.close();

    // write the solution file
    FILE* pFile;
    pFile = fopen(solutionfile,"wb");



    int nVars = values.size() / nPoints;
    const int nRestart_Vars = 5;
    int Restart_Vars[nRestart_Vars] = {-1};
    int nFields = dimension + nVars;
    Restart_Vars[0] = 535532 ; // SU2 magic numbers
    Restart_Vars[1] = nFields;
    Restart_Vars[2] = nPoints;
    Restart_Vars[3] = 1; // 1 metadata int
    Restart_Vars[4] = 8; // 8 metadata doubles


    ret = fwrite(Restart_Vars, sizeof(int), nRestart_Vars, pFile);
    CHECK( ret == nRestart_Vars);

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
        fwrite(&points[dimension*i], sizeof(double), dimension, pFile);
        fwrite(&values[i*nVars], sizeof(double),nVars,pFile);
    }

    MetaData md;
    if(md.unpack())
    {
        DPRINT(md.ext_iter);
        ret = fwrite(&md.ext_iter, sizeof(int), 1, pFile);
        CHECK( ret == 1);
        ret = fwrite(&md.aoa, sizeof(double), 1, pFile);
        CHECK( ret == 1);
        ret = fwrite(&md.sangle, sizeof(double), 1, pFile);
        CHECK( ret == 1);
        ret = fwrite(&md.initial_bc_thrust, sizeof(double), 1, pFile);
        CHECK( ret == 1);
        ret = fwrite(md.other, sizeof(double)*5, 1, pFile);
        CHECK( ret == 1);
    }


    /*--- Close the file. ---*/

    fclose(pFile);

}
void write_variables_to_ascii_restart_file( const char* solutionfile, 
                                            const char* meshfile,
                                            std::vector<double>& values)
{
    using std::FILE;
    using std::fopen;
    using std::fwrite;
    
    const int CGNS_STRING_SIZE = 33;
    int dimension = 0;
    std::size_t nPoints = 0;

    // get points from mesh file
    std::fstream mesh;
    mesh.open(meshfile);

    std::string line;

    int ret;
    CHECK(mesh.is_open());
    DPRINT(meshfile);
    while(std::getline(mesh,line) && (nPoints == 0 || dimension == 0))
    {
        // Check whether the line starts by a character
        // or if it is a comment
        if(line.size() > 0 and (isalpha(line[0])) )
        {
            ret = std::strncmp(line.c_str(),"NDIME", std::min((int)line.size(),5) );
            if( ret == 0 )
            {
                std::stringstream ss(line);
                ss >> line; // NDIME=
                ss >> dimension;
            }
            ret = std::strncmp(line.c_str(),"NPOIN", std::min((int)line.size(),5) );
            if( ret == 0 )
            {
                std::stringstream ss(line);
                ss >> line; // NPOIN=
                ss >> nPoints;
            }
        }
    }
    CHECK(ret == 0);
    std::vector< double > points;
    DPRINT(nPoints);
    double x;
    points.reserve(nPoints*dimension);
    while(std::getline(mesh,line))
    {
        std::stringstream ss(line);
        for(int i = 0 ; i < dimension; i++)
        {
            ss >> x;
            points.push_back(x);
        }
    }

    mesh.close();

    // write the solution file
    std::ofstream file;
    file.open(solutionfile);

    CHECK(file.is_open());

    int nVars = values.size() / nPoints;
    int nFields = nVars + dimension;

    std::vector<std::string> varnames;
    unpack_variable_names(varnames);

    char fieldname[CGNS_STRING_SIZE];

    file << "PointID";
    for (int iVar = 0; iVar < nFields; iVar++) 
    {
        file <<"\t";
        std::strncpy(fieldname, varnames[iVar].data(), std::min((int)varnames.size(),CGNS_STRING_SIZE));
        file << fieldname;
    }

    file <<"\n";
    file << std::setprecision(15);

    for( std::size_t i = 0 ; i < nPoints ; i++)
    {
        file << i <<"\t" ;
        for(int j = 0 ; j < dimension; j++)
            file << points[i*dimension + j] << "\t";
        for(int j = 0 ; j < nVars; j++)
            file << "\t" << values[i*nVars +j];
        file <<"\n";
    }
    MetaData md;
    if(md.unpack())
    {
        file << "EXT_ITER= " << md.ext_iter << std::endl;;
        file << "AOA= " << md.aoa << std::endl;
        file << "SIDESLIP_ANGLE= " << md.sangle << std::endl;
        file << "INITIAL_BCTHRUST= " << md.initial_bc_thrust << std::endl;
    }

    file.close();

}

void write_sol_with_scalar_vars(std::string filename,int dimension, int nVars, std::size_t nPoints, std::vector<double>& values)
{
    std::fstream file;
    file.open(filename,std::ios_base::out);
    file << "MeshVersionFormatted 2\n";
    file << "\n";
    file << "Dimension " << dimension << "\n";
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

void write_solb_with_scalar_vars(std::string filename,int dimension, int nVars, std::size_t nPoints, std::vector<double>& values)
{
    using std::FILE;
    using std::fwrite;
    using std::fopen;
    using std::fclose;

    FILE* pFile;
    pFile = fopen(filename.c_str(),"wb");
    const int code = 1;
    const int version = 2;
    const int dimension_code = 3;
    DPRINT(nVars);
    DPRINT(nPoints);
    DPRINT(values.size());
    CHECK( not values.empty());
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
namespace 
{
    int64_t read_position(std::FILE* pFile, int version)
    {
        using std::fread;

        int64_t ans = -1;
        if (version == 3)
        {
            int64_t end_position;
            int64_t res;
            res = fread(&end_position, sizeof(int64_t), 1, pFile);
            ans = end_position;
            CHECK(res == 1);
        }
        else if( version == 2)
        {
            int32_t end_position;
            int64_t res;
            res = fread(&end_position, sizeof(int32_t), 1, pFile);
            ans = end_position;
            CHECK(res == 1);
        }
        else
        {   
            CHECK(version == 2 || version == 3);
        }

        return ans;
    }
}

void read_solb_with_scalar_vars(std::string filename,int& dimension, int& nVars, std::vector<double>& values)
{
    using std::FILE;
    using std::fread;
    using std::fopen;
    using std::fclose;

    FILE* pFile;
    pFile = fopen(filename.c_str(),"rb");
    int32_t code;
    int32_t version;
    int32_t dimension_code;
    int32_t dimension_value;
    CHECK(values.empty());
    int64_t res;

    // Read file header 
    {
        res = fread(&code, sizeof(int32_t), 1, pFile);
        CHECK(res == 1);
        CHECK(code == 1);

        res = fread(&version, sizeof(int32_t), 1, pFile);
        CHECK(res == 1);
        CHECK(version == 2 || version == 3);


        // end of the upcomming section
        res = fread(&dimension_code, sizeof(int32_t), 1, pFile);
        CHECK(res == 1);
        CHECK(dimension_code == 3);

        res = read_position(pFile,version);
        CHECK(res > 0);

        res = fread(&dimension_value, sizeof(int32_t), 1, pFile);
        dimension = dimension_value;
        CHECK(res == 1);
        CHECK(dimension == 2 || dimension == 3);

    }


    // Read solution header 
    {

        int32_t keyword;
        int32_t nPoints;
        int64_t res;

        res = fread(&keyword, sizeof(int32_t), 1, pFile);
        CHECK(res == 1);
        CHECK(keyword == 62) ; // SolAtVertices = 62

        // end of the upcomming section
        res = read_position(pFile,version);
        CHECK(res > 0);


        res = fread(&nPoints, sizeof(int32_t), 1, pFile);
        CHECK(res == 1);
        CHECK(nPoints > 0);

        int32_t solutions_per_node;
        res = fread(&solutions_per_node, sizeof(int32_t), 1, pFile);
        CHECK(res == 1);
        CHECK(solutions_per_node > 0);

        nVars = solutions_per_node;

        //int end_position = ftell(pFile) + (4 + nvars)*sizeof(int) + npoints*nvars*sizeof(double);
        int32_t solutions_type; 
        for ( int32_t i = 0 ; i < nVars ; i++)
        {
            res = fread(&solutions_type, sizeof(int32_t), 1, pFile);
            CHECK(res == 1);
            CHECK( solutions_type == 1);
        }

        values.resize( nVars * nPoints);


        // read values at once
        res = fread(values.data(), sizeof(double)*values.size(), 1, pFile);
        CHECK(res == 1);
    }

    // finalize file
    {
        int32_t keyword;
        res = fread(&keyword,sizeof(int32_t),1,pFile);
        CHECK(res == 1);
        CHECK( keyword == 54) ; // = GmfEmd

        res = read_position(pFile,version);
        CHECK(res == 0);
    }


    fclose(pFile);
}
