#include <stdio.h>
#include <mpi/mpi.h>
#include <stdlib.h>

#define M 4800
#define N 600

double MA[N][M + 1], V[M + 1], MAD, R;
const int MAX = 100;
const int MIN = -100;

double some_random() {
    double range = MAX - MIN;
    double div = RAND_MAX / range;
    return MIN + (rand() / div);
}

int main(int args, char **argv) {
    int size, MyP, i, j, k, m, p;
    double t0, dt, t1, t2, t3, t4, dt1, dt2, dt3, dt4;
    MPI_Status stat;
    MPI_Init(&args, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyP);

    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++)
            if ((N * MyP + i) == j)
                MA[i][j] = 2.0 * some_random();
            else
                MA[i][j] = 1.0 * some_random();

    for (j = 0; j < M; j++)
        V[j] = -(double) (j + 1) / 2. * some_random();

    for (i = 0; i < N; i++) {
        MA[i][M] = 0;
        for (j = 0; j < M; j++)
            MA[i][M] += MA[i][j] * V[j];
    }


    t0 = MPI_Wtime();
    for (p = 0; p < size; p++) {
        for (k = 0; k < N; k++) {
            if (MyP == p) {
                MAD = 1.0 / MA[k][N * p + k];
                for (j = M; j >= N * p + k; j--)
                    MA[k][j] = MA[k][j] * MAD;
                t1 = MPI_Wtime();
                for (m = p + 1; m < size; m++)
                    MPI_Send(&MA[k][0], M + 1, MPI_DOUBLE, m, 1, MPI_COMM_WORLD);
                dt1 = MPI_Wtime() - t1;
                for (i = k + 1; i < N; i++) {
                    for (j = M; j >= N * p + k; j--)
                        MA[i][j] = MA[i][j] - MA[i][N * p + k] * MA[k][j];
                }
            } else if (MyP > p) {
                t2 = MPI_Wtime();
                MPI_Recv(V, M + 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &stat);
                dt2 = MPI_Wtime() - t2;
                for (i = 0; i < N; i++) {
                    for (j = M; j >= N * p + k; j--)
                        MA[i][j] = MA[i][j] - MA[i][N * p + k] * V[j];
                }
            }
        }
    }
    for (p = size - 1; p >= 0; p--) {
        for (k = N - 1; k >= 0; k--) {
            if (MyP == p) {
                t3 = MPI_Wtime();
                for (m = p - 1; m >= 0; m--)
                    MPI_Send(&MA[k][M], 1, MPI_DOUBLE, m, 1, MPI_COMM_WORLD);
                dt3 = MPI_Wtime() - t3;
                for (i = k - 1; i >= 0; i--)
                    MA[i][M] -= MA[k][M] * MA[i][N * p + k];
            } else {
                if (MyP < p) {
                    t4 = MPI_Wtime();
                    MPI_Recv(&R, 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &stat);
                    dt4 = MPI_Wtime() - t4;
                    for (i = N - 1; i >= 0; i--)
                        MA[i][M] -= R * MA[i][N * p + k];
                }
            }
        }
    }
    dt = MPI_Wtime() - t0;
    t0 = dt1 + dt2 + dt3 + dt4;
    printf("MyP = %d Time = %13.4e los time=%13.4e\n", MyP, dt, t0);
    MPI_Finalize();
    return 0;
}