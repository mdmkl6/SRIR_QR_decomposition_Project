#include "mpi.h"
#include <iostream>
#include <string.h>
#include <math.h>
#include "./utils.h"

#define MAX_PROCESSES 200

int main(int argc, char *argv[])
{
    int numprocs, myid, N;
    matrix *A, *Q, *mat, *p, *matTmp;
    double *vec, coef;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int displs[MAX_PROCESSES], send_counts[MAX_PROCESSES], displs2[MAX_PROCESSES], send_counts2[MAX_PROCESSES];

    // Sprawdzenie, czy podano nazwę pliku.
    if (myid == 0)
    {
        if (argc < 1)
        {
            std::cout << "Plese provide filename." << std::endl;
            MPI_Finalize();
            return 0;
        }

        // Wczytanie macierzy z pliku i sprawdzenie jej wymiarów.
        A = matrix::loadfromfile(argv[1]);
        if (A->_numcols != A->_numrows)
        {
            std::cout << "Matrix sizes not match" << std::endl;
            MPI_Finalize();
            return 0;
        }

        // Rozgłoszenie macierzy A do wszystkich procesów.
        A->broadcast_to_all();

        std::cout << "Loaded Matrix A:" << std::endl;
        A->print();
    }
    else
    {
        // Odebranie macierzy A od procesu o id=0.
        A = matrix::recive_broadcast();
    }

    N = A->_numcols;

    // Stworzenie macierzy Q jako macierzy jednostkowej o wymiarze NxN.
    Q = matrix::create_identity_matrix(N);

    int tmpN, tmpN2;

    for (int i = 0; i < N; i++)
    {
        // Obliczenie wymiarów podmacierzy, którą będzie przetwarzał każdy proces.
        tmpN = (N - i) / numprocs;
        if (myid == (numprocs - 1) && numprocs > 1)
            tmpN += (N - i) % numprocs;

        tmpN2 = (N) / numprocs;
        if (myid == (numprocs - 1) && numprocs > 1)
            tmpN2 += (N) % numprocs;

        // Usunięcie poprzednio przetwarzanej podmacierzy.
        if (i > 0 && i < N - numprocs)
        {
            delete mat;
        }

        // Stworzenie nowej podmacierzy o wymiarach tmpN x (N-i).
        mat = new matrix(tmpN, N - i);

        // Synchronizacja macierzy A i Q między procesami.
        A->synchronize_all();
        Q->synchronize_all();

        int piece = (N - i) / numprocs;
        int row = N - i;

        send_counts[0] = piece * row;
        displs[0] = 0;

        if (numprocs > 1)
        {
            for (int k = 1; k < numprocs - 1; k++)
            {
                send_counts[k] = piece * row;
                displs[k] = displs[k - 1] + piece * row;
            }
            displs[numprocs - 1] = displs[numprocs - 2] + piece * row;
            send_counts[numprocs - 1] = (((N - i) * row) - displs[numprocs - 1]);
        }

        piece = N / numprocs;
        row = N - i;
        send_counts2[0] = piece * row;
        displs2[0] = 0;
        if (numprocs > 1)
        {
            for (int k = 1; k < numprocs - 1; k++)
            {
                send_counts2[k] = piece * row;
                displs2[k] = displs2[k - 1] + piece * row;
            }
            displs2[numprocs - 1] = displs2[numprocs - 2] + piece * row;
            send_counts2[numprocs - 1] = ((N * row) - displs2[numprocs - 1]);
        }

        double x = 0;
        if (i > 0 && i < N - numprocs)
        {
            delete[] vec;
            vec = NULL;
        }
        vec = new double[N - i];
        if (myid == 0)
        {
            for (int j = 0; j < N - i; j++)
            {
                vec[j] = -A->_data[j + i][i];
                x += vec[j] * vec[j];
            }

            x = sqrt(x);

            if (vec[0] > 0)
                x = -x;
            vec[0] = vec[0] + x;
            x = 0;
            for (int j = 0; j < N - i; j++)
            {
                x += vec[j] * vec[j];
            }
            x = sqrt(x);
        }
        MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(vec, N - i, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (x > 0)
        {
            if (myid == 0)
            {
                for (int j = 0; j < N - i; j++)
                {
                    vec[j] /= x;
                }
            }

            MPI_Bcast(vec, N - i, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (i > 0 && i < N - numprocs)
            {
                delete p;
            }
            p = new matrix(N - i, N - i);

            p->scatter(send_counts, displs, mat, tmpN * (N - i));

            for (int k = 0; k < send_counts[myid] / (N - i); k++)
            {
                for (int l = 0; l < N - i; l++)
                {
                    if ((k + (displs[myid] / (N - i))) == l)
                        mat->_data[k][l] = 1 - 2 * vec[k + displs[myid] / (N - i)] * vec[l];
                    else
                        mat->_data[k][l] = -2 * vec[k + displs[myid] / (N - i)] * vec[l];
                }
            }
            mat->gather(send_counts[myid], p, send_counts, displs);
            p->synchronize_all();

            // R
            for (int k = 0; k < send_counts[myid] / (N - i); k++)
            {
                for (int l = 0; l < N - i; l++)
                {
                    double tm = 0;
                    for (int m = i; m < N; m++)
                    {
                        tm += p->_data[k + displs[myid] / (N - i)][m - i] * A->_data[m][l + i];
                    }
                    mat->_data[k][l] = tm;
                }
            }

            if (i > 0 && i < N - numprocs)
            {
                delete matTmp;
            }
            matTmp = new matrix(N - i, N - i);
            mat->gather(send_counts[myid], matTmp, send_counts, displs);

            if (myid == 0)
            {
                for (int k = i; k < N; k++)
                {
                    for (int l = i; l < N; l++)
                    {
                        A->_data[k][l] = matTmp->_data[k - i][l - i];
                    }
                }
            }

            // Q
            if (i > 0 && i < N - numprocs)
            {
                delete mat;
            }
            mat = new matrix(tmpN2, N - i);

            for (int k = 0; k < send_counts2[myid] / (N - i); k++)
            {
                for (int l = 0; l < N - i; l++)
                {
                    double tm = 0;
                    for (int m = i; m < N; m++)
                    {
                        tm += Q->_data[k + displs2[myid] / (N - i)][m] * p->_data[m - i][l];
                    }
                    mat->_data[k][l] = tm;
                }
            }

            if (i > 0 && i < N - numprocs)
            {
                delete matTmp;
            }
            matTmp = new matrix(N, N - i);
            mat->gather(send_counts2[myid], matTmp, send_counts2, displs2);

            if (myid == 0)
            {
                for (int k = 0; k < N; k++)
                {
                    for (int l = i; l < N; l++)
                    {
                        Q->_data[k][l] = matTmp->_data[k][l - i];
                    }
                }
            }
        }
    }

if (myid == 0)
{
    // Jeśli identyfikator bieżącego procesu to 0, to wyświetl rozwiązanie.
    std::cout << "Solution" << std::endl;

    // Wyświetl macierz Q.
    std::cout << "Matrix Q:" << std::endl;
    Q->print();

    // Wyświetl macierz R.
    std::cout << "Matrix R:" << std::endl;
    A->print();

    if (argc >= 2)
        {
            std::ofstream outfile(argv[2]);
            if (outfile.is_open())
            {
                Q->save_to_file(outfile);
                outfile<<std::endl;
                A->save_to_file(outfile);
                outfile.close();
            }
        }
}

// Zwolnij pamięć, która została dynamicznie przydzielona dla obiektów A, Q, p, matTmp i vec.
delete A;
delete Q; 
delete p;
delete matTmp;
delete[] vec;

// Zakończ pracę biblioteki MPI.
MPI_Finalize();

return 0;
}