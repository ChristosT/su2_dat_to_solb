#ifndef SOLUTION_IO_H
#define SOLUTION_IO_H
#include <vector>
#include <string>
/**
 * @brief read all variable values on every point from a su2 binary .dat file
 *
 * @param filename name of the solution file
 * @param[out] values  flat table of variable values in the form [ var0 for p0, var1, for p0,...]
 * @param[out] varnames a vector of names of variables appearing in the file. index matches with column ids
 * @param[out] nPoints number of points in the solution
 * @param[out] dimension the dimension of the points (2 or 3)
 */
void read_variables_from_binary_restart_file( const char* filename, 
                                              std::vector<double>& values,
                                              std::vector<std::string>& varnames,
                                              std::size_t& nPoints,
                                              int& dimension);



/**
 * @brief write a binary/ascii restart file for SU2
 *
 * @param solutionfile name of the file to be created
 * @param meshfile name of the matching SU2 mesh file
 * @param values  a flat table of size nvars*nPoints of variable values in the form [ var0 for p0, var1, for p0,...]
 */
void write_variables_to_binary_restart_file( const char* solutionfile, 
                                             const char* meshfile,
                                             std::vector<double>& values);

void write_variables_to_ascii_restart_file( const char* solutionfile, 
                                             const char* meshfile,
                                             std::vector<double>& values);

/**
 * @brief write an ascii(sol)/binary(solb) Gamma Mesh Format (GMF) solution file of \p nvar scalar variables per vertex
 *
 * @param basefilename name of the file wqith no suffix
 * @param dimension the dimension of the points (2 or 3)
 * @param nVars number of variables(per point) contained in the values vector
 * @param nPoints number of points in the solution
 * @param values  flat table of variable values in the form [ var0 for p0, var1, for p0,...]
 *
 * @note uses SolAtVertices keyword
 */
void write_sol_with_scalar_vars(std::string basefilename,int dimension, int nVars, std::size_t nPoints, std::vector<double>& values);
void write_solb_with_scalar_vars(std::string basefilename,int dimension, int nVars,std::size_t nPoints,  std::vector<double>& values);

/**
 * @brief read a binary(solb) Gamma Mesh Format (GMF) solution file of \p nvar scalar variables per vertex
 *
 * @param filename name of the file solution file
 * @param dimension the dimension of the points (2 or 3)
 * @param nvars number of variables(per point) contained in the values vector
 * @param values  flat table of variable values in the form [ var0 for p0, var1, for p0,...]
 *
 * @note assumes  SolAtVertices is used
 */
void read_solb_with_scalar_vars(std::string filename,int& dimension, int& nvars, std::vector<double>& values);
#endif /* SOLUTION_IO_H */
