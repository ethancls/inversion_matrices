#include "mat.h"

int main(int argc, char **argv)
{

    if (argc != 3)
    {

        printf("Error, 2 matrix needed\n");
        exit(EXIT_FAILURE);
    }

    uint32_t tailleA, tailleB;
    double **A = readFromFile(argv[1], &tailleA);
    double **B = readFromFile(argv[2], &tailleB);

    if (tailleA == tailleB)
    {
        double **C = multiply(A, B, tailleA);
        printMatrix(C, tailleA);
        freeMatrix(C, tailleA);
    }
    else
    {
        printf("Error with matrix size\n");
        exit(EXIT_FAILURE);
    }

    freeMatrix(A, tailleA);
    freeMatrix(B, tailleB);

    return EXIT_SUCCESS;
}