#include <iostream> 
#include <fstream>
#include <Eigen/Cholesky>
#include <unsupported/Eigen/MPRealSupport>

using namespace Eigen;

//using MatrixXmp = Matrix<mpfr::mpreal, Dynamic, Dynamic>;
//using VectorXmp = Matrix<mpfr::mpreal, Dynamic, 1>;
//using namespace mpfr;
typedef Matrix<mpfr::mpreal, Dynamic, Dynamic> MatrixXmp;

bool check_matrix_psd(MatrixXmp & A)
{
    //check if the matrix is positive definite
    LLT<MatrixXmp> lltOfA(A);
    if (lltOfA.info() == NumericalIssue)
    {
        std::cerr << "Error: matrix is not positive definite" << std::endl;
        return false;
    }
    std::cout << "Matrix is positive definite" << std::endl;

    return true;
}

void read_matrices_from_file_and_check(std::ifstream & file)
{
    int linenum = 1;
    // while the file has more lines, read each line of the file
    while (file)
    {
       // read a line from the file
        std::cout << "Reading line " << linenum++ << std::endl;
        std::string line;
        std::getline(file, line);
        // if line is empty or all spaces, skip it
        if (line.empty() || line.find_first_not_of(" \t\n\v\f\r") == std::string::npos)
        {
            continue;
        }
        //the line is a matrix written out in row-major order
        //parse the line to get the matrix size 
        // the line starts with a character '[' 
        assert(line[0] == '[');
        int pos = 1;
        std::vector< mpfr::mpreal > matrix_data;
        // split the line into tokens with separators that include ' ' and ';' 
        while (pos < line.size())
        {
            // find the next token
            int next_pos = line.find_first_of(" ;]", pos);
            // extract the token
            if (next_pos >= line.size())
            {
                break;
            }

            if (next_pos == pos){
                pos++;
                continue;
            }
            std::string token = line.substr(pos, next_pos - pos);
            // convert the token to a number
            mpfr::mpreal number(token);
//            std::cout << "n: " << number << std::endl;
            // assert that the numbe is a valid number
            assert (mpfr::isfinite(number));
            // add the number to the matrix data
            matrix_data.push_back(number);
            // move to the next position
            pos = next_pos + 1;
        }

        if (matrix_data.size() == 1)
        {
            // the matrix is a scalar
            std::cout << "Matrix is a scalar" << std::endl;
            if (matrix_data[0] < 0)
            {
                std::cerr << "Error: matrix is not positive definite" << std::endl;
            }
            continue;
        }
        int n = sqrt(matrix_data.size());

        // check if the matrix data is square
        std::cout << "\t Matrix size is: " << n << " I have loaded " << matrix_data.size()  << " numbers" << std::endl;
       
        assert(n * n == matrix_data.size());
        // create a matrix from the data
        MatrixXmp A(n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A(i, j) = matrix_data[i * n + j];
            }
        }
        // print the matrix out 
        //std::cout << "Matrix is: " << std::endl;
        //std::cout << A << std::endl;

        //check if the matrix is positive definite
        check_matrix_psd(A);
    }
}

int main(int argc, char *argv[])
{
    mpfr::mpreal::set_default_prec(512);
    // Read in the matrix from a file 
    //open a file for reading 
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    std::ifstream file(argv[1]);
    //check if the file is open
    if (!file.is_open())
    {
        std::cerr << "Error: file not found" << std::endl;
        return 1;
    }
    //print("opening file: ", argv[1]);
    std::cout << "opening file: " << argv[1] << std::endl;
    //read the matrix size
    read_matrices_from_file_and_check(file);
    //close the file
    file.close();
   
   
    return 0;
}