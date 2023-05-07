#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Funkcja do wczytywania macierzy z pliku
std::vector<std::vector<double>> read_matrix(std::ifstream &infile)
{
    int rowlen = 0;
    std::vector<std::vector<double>> matrix;
    std::string line;
    while (std::getline(infile, line))
    {
        std::vector<double> row;
        std::istringstream iss(line);
        double value;
        while (iss >> value)
        {
            row.push_back(value);
        }
        if (rowlen == 0)
            rowlen = row.size();
        matrix.push_back(row);
        if (matrix.size() >= rowlen)
        {
            std::getline(infile, line);
            break;
        }
    }
    return matrix;
}

// Funkcja do wymnażania macierzy
std::vector<std::vector<double>> matrix_multiply(const std::vector<std::vector<double>> &A,
                                                 const std::vector<std::vector<double>> &B)
{
    const size_t n = A.size();
    const size_t m = B[0].size();
    const size_t p = B.size();
    std::vector<std::vector<double>> C(n, std::vector<double>(m, 0.0));
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            for (size_t k = 0; k < p; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double meanError(const std::vector<std::vector<double>> &matrix1, const std::vector<std::vector<double>> &matrix2)
{
    double sumError = 0.0;
    const size_t m = matrix1.size();
    const size_t n = matrix1[0].size();
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            sumError += pow(matrix1[i][j] - matrix2[i][j], 2.0);
        }
    }
    const double meanError = sumError / static_cast<double>(m * n);
    return meanError;
}

void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& elem : row) {
            std::cout << std::setw(9) << std::setfill(' ') << std::setprecision(5) << std::fixed << elem << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Not enaugh arguments!" << std::endl;
        return 1;
    }

    std::ifstream infile(argv[1]);
    if (!infile.is_open())
    {
        std::cout << "Failed to open file: " << argv[1] << std::endl;
        return 1;
    }
    std::vector<std::vector<double>> A = read_matrix(infile);
    infile.close();

    infile = std::ifstream(argv[2]);
    if (!infile.is_open())
    {
        std::cout << "Failed to open file: " << argv[2] << std::endl;
        return 1;
    }
    std::vector<std::vector<double>> Q = read_matrix(infile);
    std::vector<std::vector<double>> R = read_matrix(infile);
    infile.close();

    // Wymnożenie macierzy Q i R, aby otrzymać macierz QR
    std::vector<std::vector<double>> QR = matrix_multiply(Q, R);

    // Obliczenie średniego błędu między macierzą A a macierzą QR
    double error = meanError(A, QR);
    std::cout << "Mean error: " << error << std::endl;

    std::cout << "Q*R matrix:" << std::endl;
    printMatrix(QR);

    return 0;
}