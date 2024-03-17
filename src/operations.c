#include "mat.h"

// Allocate a 2D array for double matrix
double **dataMatrix(uint32_t n)
{
    double **mat = (double **)calloc(n, sizeof(double *));
    for (uint32_t i = 0; i < n; i++)
    {
        mat[i] = (double *)calloc(n, sizeof(double));
    }
    return mat;
}

// Free a matrix (2D array)
void freeMatrix(double **mat, uint32_t n)
{
    for (uint32_t i = 0; i < n; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

// Print matrix on terminal
void printMatrix(double **A, uint32_t n)
{
    for (uint32_t i = 0; i < n; i++)
    {
        printf("| ");
        for (uint32_t j = 0; j < n; j++)
        {
            printf("%.4lf  ", A[i][j]);
        }
        printf("|\n");
    }
}

// Print matrix to a file
void printMatrixToFile(const char *f, double **A, uint32_t n)
{
    FILE *file = fopen(f, "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file: %s\n", f);
        return;
    }

    for (uint32_t i = 0; i < n; i++)
    {
        fprintf(file, "| ");
        for (uint32_t j = 0; j < n; j++)
        {
            fprintf(file, "%.5lf  ", A[i][j]);
        }
        fprintf(file, "|\n");
    }

    fclose(file);
}

// Function to print the matrix in the specified format
void printMatrixFormatted(double **A, uint32_t n)
{
    printf("[");
    for (uint32_t i = 0; i < 3; i++)
    {
        printf("[");
        for (uint32_t j = 0; j < 3; j++)
        {
            printf("%.8lf", A[i][j]);
            if (j < 2)
                printf(" ");
        }
        printf("...");
        for (uint32_t j = n - 3; j < n; j++)
        {
            printf("%.8lf", A[i][j]);
            if (j < n - 1)
                printf(" ");
        }
        printf("]");
        if (i < 2)
            printf("\n");
    }
    printf("\n...\n ");
    for (uint32_t i = n - 3; i < n; i++)
    {
        printf("[");
        for (uint32_t j = 0; j < 3; j++)
        {
            printf("%.8lf", A[i][j]);
            if (j < 3)
                printf(" ");
        }
        printf("...");
        for (uint32_t j = n - 3; j < n; j++)
        {
            printf("%.8lf", A[i][j]);
            if (j < n - 1)
                printf(" ");
        }
        printf("]");
        if (i < n - 1)
            printf("\n");
    }
    printf("]\n");
}

// Add matrix A to B
double **add(double **A, double **B, uint32_t n)
{
    double **C = dataMatrix(n);
    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// Subtract matrix B from A
double **subtract(double **A, double **B, uint32_t n)
{
    double **C = dataMatrix(n);
    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

// Return the negative matrix (same as -1 scalar multiplication)
double **negate(double **A, uint32_t n)
{
    double **C = dataMatrix(n);
    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            C[i][j] = -A[i][j];
        }
    }
    return C;
}

// Transpose matrix A to tA
double **transpose(double **A, uint32_t n)
{
    double **B = dataMatrix(n);
    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            B[j][i] = A[i][j];
        }
    }
    return B;
}

// Standard multiplication
double **multiply(double **A, double **B, uint32_t n)
{
    double **C = dataMatrix(n);
    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            C[i][j] = 0.0;
            for (uint32_t k = 0; k < n; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Read a matrixfrom a file
double **readFromFile(const char *fileName, uint32_t *n)
{
    FILE *file = fopen(fileName, "r");

    if (file == NULL)
    {
        perror("Error reading file");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%u", n);

    double **matrix = dataMatrix(*n);

    for (uint32_t i = 0; i < *n; i++)
    {
        for (uint32_t j = 0; j < *n; j++)
        {
            fscanf(file, "%lf", &matrix[i][j]);
        }
    }

    fclose(file);
    return matrix;
}

// Combine quadrants into one matrix of 'size'
void combineMatrices(double **result, double **C11, double **C12, double **C21, double **C22, uint32_t size)
{
    for (uint32_t i = 0; i < size / 2; i++)
    {
        for (uint32_t j = 0; j < size / 2; j++)
        {
            result[i][j] = C11[i][j];
            result[i][j + size / 2] = C12[i][j];
            result[i + size / 2][j] = C21[i][j];
            result[i + size / 2][j + size / 2] = C22[i][j];
        }
    }
}