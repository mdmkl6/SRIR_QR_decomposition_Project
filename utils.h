#include "mpi.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

class matrix
{
public:
    double **_data;
    int _numrows, _numcols;

    matrix(int numrows, int numcols)
    {
        _numrows = numrows;
        _numcols = numcols;
        double *buffor = new double[numrows * numcols];
        _data = new double *[numrows];
        for (int i = 0; i < numrows; ++i)
            _data[i] = buffor + i * numcols;
    }

    matrix(int numrows, int numcols, double *buffor)
    {
        _numrows = numrows;
        _numcols = numcols;
        _data = new double *[numrows];
        for (int i = 0; i < numrows; ++i)
            _data[i] = buffor + i * numcols;
    }

    ~matrix()
    {
        delete[] _data[0];
        _numrows = 0;
        _numcols = 0;
    }

    void scatter(int *sendcounts, int *displs, matrix *target, int recvcount)
    {
        MPI_Scatterv(&_data[0][0], sendcounts, displs, MPI_DOUBLE, target->_data[0], recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void gather(int sendcounts, matrix *target, int *recvcount, int *displs)
    {
        MPI_Gatherv(&_data[0][0], sendcounts, MPI_DOUBLE, target->_data[0], recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void print()
    {
        for (int i = 0; i < _numrows; i++)
        {
            for (int j = 0; j < _numcols; j++)
            {
                std::cout <<  std::setw(9) << std::setfill(' ') << std::setprecision(5) << std::fixed << _data[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void broadcast_to_all()
    {
        MPI_Bcast(&_numrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&_numcols, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(_data[0], _numrows * _numcols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void synchronize_all()
    {
        MPI_Bcast(_data[0], _numrows * _numcols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    static matrix *loadfromfile(std::string filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Failed to open file " << filename << std::endl;
            MPI_Finalize();
            exit(1);
        }

        std::vector<std::vector<double>> matrix_data;

        std::string line;
        while (std::getline(file, line))
        {
            std::vector<double> row_data;
            std::stringstream line_stream(line);
            double val;
            while (line_stream >> val)
            {
                row_data.push_back(val);
            }
            matrix_data.push_back(row_data);
        }
        file.close();

        int numrows = matrix_data.size();
        int numcols = matrix_data[0].size();
        matrix *mat = new matrix(numrows, numcols);

        for (int i = 0; i < numrows; ++i)
        {
            for (int j = 0; j < numcols; ++j)
            {
                mat->_data[i][j] = matrix_data[i][j];
            }
        }

        return mat;
    }

    static matrix *recive_broadcast()
    {
        int numrows, numcols;
        MPI_Bcast(&numrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numcols, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double *buffor = new double[numrows * numcols];
        MPI_Bcast(buffor, numrows * numcols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        return new matrix(numrows, numcols, buffor);
    }

    static matrix *create_identity_matrix(int n)
    {
        matrix *mat = new matrix(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    mat->_data[i][j] = 1;
                else
                    mat->_data[i][j] = 0;
            }
        return mat;
    }
};