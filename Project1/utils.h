#include "mpi.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

class matrix
{
public:
    double **_data; // wskaźnik na dwuwymiarową tablicę double'ów
    int _numrows, _numcols; // liczba wierszy i liczba kolumn macierzy

    // Konstruktor
    matrix(int numrows, int numcols)
    {
        _numrows = numrows;
        _numcols = numcols;
        double *buffor = new double[numrows * numcols]; // bufor na dane macierzy
        _data = new double *[numrows]; // alokacja pamięci na wskaźniki do wierszy
        for (int i = 0; i < numrows; ++i)
            _data[i] = buffor + i * numcols; // wskaźniki wierszy ustawione na odpowiednie miejsca w buforze
    }
    // Konstruktor z danymi początkowymi
    matrix(int numrows, int numcols, double *buffor)
    {
        _numrows = numrows;
        _numcols = numcols;
        _data = new double *[numrows];
        for (int i = 0; i < numrows; ++i)
            _data[i] = buffor + i * numcols;
    }

    // Destruktor
    ~matrix()
    {
        delete[] _data[0]; // zwolnienie bufora na dane
        _numrows = 0;
        _numcols = 0;
    }

    // Metoda do rozsyłania danych macierzy pomiędzy procesami.
    // sendcounts - liczba elementów, która ma zostać wysłana do każdego z procesów
    // displs - przesunięcie początku danych w buforze sendbuf dla każdego procesu
    // target - wskaźnik do macierzy, która otrzyma rozsyłane dane
    // recvcount - liczba elementów, które zostaną odebrane przez proces
    void scatter(int *sendcounts, int *displs, matrix *target, int recvcount)
    {
        MPI_Scatterv(&_data[0][0], sendcounts, displs, MPI_DOUBLE, target->_data[0], recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Metoda do zbierania danych macierzy od procesów.
    // sendcounts - liczba elementów, która ma zostać wysłana przez każdy z procesów
    // target - wskaźnik do macierzy, która otrzyma zebrane dane
    // recvcount - liczba elementów, które zostaną odebrane przez proces
    // displs - przesunięcie początku danych w buforze recvbuf dla każdego procesu
    void gather(int sendcounts, matrix *target, int *recvcount, int *displs)
    {
        MPI_Gatherv(&_data[0][0], sendcounts, MPI_DOUBLE, target->_data[0], recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Metoda do wyświetlania zawartości macierzy na ekranie.
    void print()
    {
        // iterujemy po wierszach macierzy
        for (int i = 0; i < _numrows; i++)
        {
            // iterujemy po kolumnach macierzy
            for (int j = 0; j < _numcols; j++)
            {
                // wyświetlamy każdy element z dokładnością do 5 miejsc po przecinku i odstępem 9 znaków
                std::cout << std::setw(9) << std::setfill(' ') << std::setprecision(5) << std::fixed << _data[i][j] << " ";
            }
            // przechodzimy do następnego wiersza
            std::cout << std::endl;
        }
        // wstawiamy dodatkową linię dla czytelności
        std::cout << std::endl;
    }
    void save_to_file(std::ofstream &file)
    {
        for (int i = 0; i < _numrows; i++)
        {
            for (int j = 0; j < _numcols; j++)
            {
                file << _data[i][j] << " ";
            }
            file << std::endl;
        }
    }
    void broadcast_to_all()
    {
        // broadcastujemy liczbę wierszy
        MPI_Bcast(&_numrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // broadcastujemy liczbę kolumn
        MPI_Bcast(&_numcols, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // broadcastujemy całą macierz
        MPI_Bcast(_data[0], _numrows * _numcols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Metoda do synchronizacji macierzy między procesami.
    void synchronize_all()
    {
        // broadcastujemy całą macierz
        MPI_Bcast(_data[0], _numrows * _numcols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Metoda do wczytywania macierzy z pliku.
    static matrix *loadfromfile(std::string filename)
    {
        // otwieramy plik
        std::ifstream file(filename);
        // sprawdzamy, czy udało się otworzyć plik
        if (!file.is_open())
        {
            std::cerr << "Failed to open file " << filename << std::endl;
            MPI_Finalize();
            exit(1);
        }

        // deklarujemy wektor dwuwymiarowy przechowujący elementy macierzy
        std::vector<std::vector<double>> matrix_data;

        // wczytujemy plik linia po linii
        std::string line;
        while (std::getline(file, line))
        {
            // deklarujemy wektor przechowujący wartości w danym wierszu
            std::vector<double> row_data;
            // dzielimy linię na wartości oddzielone spacją
            std::stringstream line_stream(line);
            double val;
            while (line_stream >> val)
            {
                row_data.push_back(val);
            }
            // dodajemy wiersz do macierzy
            matrix_data.push_back(row_data);
        }
        // zamykamy plik
        file.close();

        // tworzymy nową macierz na podstawie wczytanych wartości
        int numrows = matrix_data.size(); // liczba wierszy macierzy
        int numcols = matrix_data[0].size(); // liczba kolumn macierzy
        matrix *mat = new matrix(numrows, numcols); // utworzenie obiekt

        for (int i = 0; i < numrows; ++i)
        {
            for (int j = 0; j < numcols; ++j)
            {
                mat->_data[i][j] = matrix_data[i][j];
            }
        }

        return mat;
    }
    
    // Metoda, która odbiera i rozgłasza macierz przy użyciu MPI
    // Zwraca wskaźnik na nowo utworzoną macierz
    static matrix *recive_broadcast()
    {
        int numrows, numcols;
        
        // Odbieranie liczby wierszy i kolumn macierzy
        MPI_Bcast(&numrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numcols, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Tworzenie bufora na dane macierzy
        double *buffor = new double[numrows * numcols];

        // Odbieranie danych macierzy
        MPI_Bcast(buffor, numrows * numcols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Tworzenie nowej macierzy na podstawie otrzymanych danych i zwracanie wskaźnika na nią
        return new matrix(numrows, numcols, buffor);
    }

    // Metoda, która tworzy macierz jednostkową o wymiarach nxn
    // Zwraca wskaźnik na nowo utworzoną macierz
    static matrix *create_identity_matrix(int n)
    {
        // Tworzenie nowej macierzy o wymiarach nxn
        matrix *mat = new matrix(n, n);

        // Wypełnianie macierzy wartościami
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    mat->_data[i][j] = 1; // wartość 1 na przekątnej
                else
                    mat->_data[i][j] = 0; // wartość 0 poza przekątną
            }

        // Zwracanie wskaźnika na nową macierz jednostkową
        return mat;
    }
};