#include <upcxx/upcxx.hpp>
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
    void scatter(int myid, matrix *target, int recvstart, int count)
    {
        upcxx::global_ptr<double> matrix;
        if(myid==0)
        {
            matrix = upcxx::new_array<double>(_numrows * _numcols);
            double *local_matrix = matrix.local();
            for (int i = 0; i < _numrows * _numcols; i++) {
                local_matrix[i] = _data[0][i];
            }
            matrix = upcxx::broadcast(matrix, 0).wait();
            for(int i = 0; i < count; i++)
                target->_data[0][i]=_data[0][recvstart+i];
        }
        else{
            matrix = upcxx::broadcast(matrix, 0).wait();
            for(int i = 0; i < count; i++)
                target->_data[0][i]=upcxx::rget(matrix+recvstart+i).wait();
        }
        //MPI_Scatterv(&_data[0][0], sendcounts, displs, MPI_DOUBLE, target->_data[0], recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Metoda do zbierania danych macierzy od procesów.
    // sendcounts - liczba elementów, która ma zostać wysłana przez każdy z procesów
    // target - wskaźnik do macierzy, która otrzyma zebrane dane
    // recvcount - liczba elementów, które zostaną odebrane przez proces
    // displs - przesunięcie początku danych w buforze recvbuf dla każdego procesu
    void gather(int myid, matrix *target, int sendstart, int count)
    {
        std::cout<<sendstart<<std::endl;
        std::cout<<count<<std::endl;
        upcxx::global_ptr<double> matrix;
        if(myid==0)
        {
            matrix = upcxx::new_array<double>(_numrows * _numcols);
            matrix = upcxx::broadcast(matrix, 0).wait();

            int i = 0;
            for(; i < count; i++)
                target->_data[0][sendstart+i]=_data[0][i];

            upcxx::barrier();
            double *local_matrix = matrix.local();
            for(; i < _numrows * _numcols; i++)
                target->_data[0][sendstart+i]=local_matrix[i];
        }
        else{
            matrix = upcxx::broadcast(matrix, 0).wait();
            for(int i = 0; i < count; i++)
                upcxx::rput(_data[0][i], matrix+sendstart+i).wait();
            upcxx::barrier();
        }
        
        //MPI_Gatherv(&_data[0][0], sendcounts, MPI_DOUBLE, target->_data[0], recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
        upcxx::broadcast(&_numrows, 1, 0).wait();
        // broadcastujemy liczbę kolumn
        upcxx::broadcast(&_numcols, 1, 0).wait();
        // broadcastujemy całą macierz
        upcxx::broadcast(_data[0], _numrows * _numcols, 0).wait();
    }

    // Metoda do synchronizacji macierzy między procesami.
    void synchronize_all()
    {
        // broadcastujemy całą macierz
        upcxx::broadcast(_data[0], _numrows * _numcols, 0).wait();
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
            upcxx::finalize();
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
        upcxx::broadcast(&numrows, 1, 0).wait();
        upcxx::broadcast(&numcols, 1, 0).wait();

        // Tworzenie bufora na dane macierzy
        double *buffor = new double[numrows * numcols];

        // Odbieranie danych macierzy
        upcxx::broadcast(buffor, numrows * numcols, 0).wait();

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