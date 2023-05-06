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

    if (myid == 0)
    {
        if (argc < 1)
        {
            std::cout << "Plese provide filename." << std::endl;
            MPI_Finalize();
            return 0;
        }
        A = matrix::loadfromfile(argv[1]);
        if (A->_numcols != A->_numrows)
        {
            std::cout << "Matrix sizes not match" << std::endl;
            MPI_Finalize();
            return 0;
        }
        A->broadcast_to_all();

        std::cout << "Loaded Matrix A:" << std::endl;
        A->print();
    }
    else
    {
        A = matrix::recive_broadcast();
    }
    N = A->_numcols;
    Q = matrix::create_identity_matrix(N);

    int tmpN, tmpN2;

    for (int i = 0; i < N; i++)
    {
        tmpN = (N - i) / numprocs;
        if (myid == (numprocs - 1) && numprocs > 1)
            tmpN += (N - i) % numprocs;

        tmpN2 = (N) / numprocs;
        if (myid == (numprocs - 1) && numprocs > 1)
            tmpN2 += (N) % numprocs;

        if (i > 0 && i < N - numprocs)
        {
            delete mat;
        }
        mat = new matrix(tmpN, N - i);

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

        float x = 0;
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
        std::cout << "Solution" << std::endl;
        std::cout << "Matrix Q:" << std::endl;
        Q->print();
        std::cout << "Matrix R:" << std::endl;
        A->print();
    }

    delete A;
    delete Q; 
    delete p;
    delete matTmp;
    delete[] vec;
    MPI_Finalize();
    return 0;
}