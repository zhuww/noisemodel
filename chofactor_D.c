/*
 * CULA Example: basicUsage
 *
 * This example is meant to show a typical program flow when using CULA.
 *
 * 1. Allocate a matrix on the host
 * 2. Initialize CULA
 * 3. Initialize the A matrix to zero
 * 4. Call culaSgeqrf (QR factorization) on the matrix
 * 5. Call culaShutdown
 *
 * After each CULA operation, the status of CULA is checked.  On failure, an
 * error message is printed and the program exits.
 *
 * Note: this example performs the QR factorization on a matrix of zeros, the
 * result of which is a matrix of zeros, due to the omission of ones across the
 * diagonal of the upper-diagonal unitary Q.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cula_lapack.h>
#include <time.h>
#include <stdlib.h>


void checkStatus(culaStatus status)
{
    char buf[256];

    if(!status)
        return;

    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);

    culaShutdown();
    exit(EXIT_FAILURE);
}

void chofactor(int N, culaDouble* A)
{
    int i, j;
    culaStatus status;

    status = culaInitialize();
    checkStatus(status);
    status = culaDpotrf('U', N, A, N);
    checkStatus(status);
    culaShutdown();

    for (i=0;i<N;i++){
        for (j=0;j<i;j++){
            /**(Out+j*N+i) = *(A+j*N+i);*/
            *(A+j*N+i) = 0.;
        };
    }
}

void chosolve(int N, int Nrhs, culaDouble* C, int ldc, culaDouble* b, int ldb)
{
    culaStatus status;
    status = culaInitialize();
    checkStatus(status);
    status = culaDpotrs('U', N, Nrhs, C, ldc, b, ldb);
    checkStatus(status);
    culaShutdown();
}

int main(int argc, char** argv)
{
    int i, j, N, M;
    int INFO;
    culaDouble r;
    N = 10;
    M = 10;
    srand(time(NULL));

    culaDouble* A = NULL;
    printf("Allocating Matrices\n");
    A = (culaDouble*)malloc(N*N*sizeof(culaDouble));
    memset(A, 0, N*N*sizeof(culaDouble));
    for (i=0;i<N;i++)
    {
        for (j=i;j<N;j++){
            r = ((culaDouble) (rand() % 10))/30.;
            *(A+(i*N)+j) = r;
            *(A+(j*N)+i) = r;
        }
        *(A+(i*N)+i) = 1.;
    };

    /*for (i=0;i<N;i++){*/
        /*for (j=0;j<N;j++){*/
            /*printf("%0.2f ", *(A+(j*N)+i));*/
        /*};*/
        /*printf("\n");*/
    /*}*/

    chofactor(N, A);
    printf("\n");

    culaDouble* b = NULL;
    b = (culaDouble*)malloc(N*M*sizeof(culaDouble));
    memset(b, 0, N*M*sizeof(culaDouble));
    for (i=0;i<N;i++){
        for (j=0;j<M;j++){
            r = ((culaDouble) (rand() % 10))/10;
            *(b+(i*M)+j) = r;
        }
    }
    chosolve(N, N, A, N, b, M);

    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            printf("%0.2f ", *(b+(j*N)+i));
        };
        printf("\n");
    }

    /*free(A);*/

    return EXIT_SUCCESS;
}

