#include "mat.h"

// Basic Inverse with standard multiply
double **inverse(double **A, uint32_t size)
{
    // Base case : 1x1 matrix
    if (size == 1)
    {
        if (A[0][0] == 0)
        {
            printf("Not inversible !\n");
            exit(EXIT_FAILURE);
        }
        else
        {
            double **result = dataMatrix(1);
            result[0][0] = 1.0 / A[0][0];
            return result;
        }
    }
    else
    {

        // Slit matrix
        double **B = dataMatrix(size / 2);
        double **CT = dataMatrix(size / 2);
        double **C = dataMatrix(size / 2);
        double **D = dataMatrix(size / 2);

        for (uint32_t i = 0; i < size / 2; i++)
        {
            for (uint32_t j = 0; j < size / 2; j++)
            {
                B[i][j] = A[i][j];
                CT[i][j] = A[i][j + size / 2];
                C[i][j] = A[i + size / 2][j];
                D[i][j] = A[i + size / 2][j + size / 2];
            }
        }

        // Recursive inversion of B
        double **B_1 = inverse(B, size / 2);

        // Allocate memory for other matrices
        double **CB_1 = multiply(C, B_1, size / 2);
        double **B_1CT = transpose(CB_1, size / 2);
        double **CB_1CT = multiply(CB_1, CT, size / 2);
        double **S = subtract(D, CB_1CT, size / 2);

        // Recursive inversion of S
        double **S_1 = inverse(S, size / 2);

        // S−1CB−1 B−1CT S−1CB−1 | B−1 + B−1CT S−1CB−1.
        double **S_1CB_1 = multiply(S_1, CB_1, size / 2);
        double **B_1CTS_1 = transpose(S_1CB_1, size / 2);
        double **B_1CTS_1CB_1 = multiply(B_1CTS_1, CB_1, size / 2);

        double **tempB = add(B_1, B_1CTS_1CB_1, size / 2);
        double **tempCT = negate(B_1CTS_1, size / 2);
        double **tempC = negate(S_1CB_1, size / 2);

        // ATA_1
        double **ATA_1 = dataMatrix(size);
        combineMatrices(ATA_1, tempB, tempCT, tempC, S_1, size);

        freeMatrix(B, size / 2);
        freeMatrix(CT, size / 2);
        freeMatrix(C, size / 2);
        freeMatrix(D, size / 2);
        freeMatrix(B_1, size / 2);
        freeMatrix(CB_1, size / 2);
        freeMatrix(B_1CT, size / 2);
        freeMatrix(CB_1CT, size / 2);
        freeMatrix(S, size / 2);
        freeMatrix(S_1, size / 2);
        freeMatrix(S_1CB_1, size / 2);
        freeMatrix(B_1CTS_1, size / 2);
        freeMatrix(B_1CTS_1CB_1, size / 2);
        freeMatrix(tempB, size / 2);
        freeMatrix(tempCT, size / 2);
        freeMatrix(tempC, size / 2);

        return ATA_1;
    }
}

int main(int argc, char *argv[])
{
    srand(time(NULL));

    uint32_t n;
    clock_t start, end;
    double cpu_time_used;

    printf("1. Use existing matrix file\n");
    printf("2. Generate a matrix\n");
    printf("Enter your choice (1 or 2): ");

    int choice, generate;
    double **A;
    char filename[100];
    scanf("%d", &choice);

    switch (choice)
    {
    case 1:
        printf("\nEnter the name of the matrix file: ");
        scanf("%s", filename);
        A = readFromFile(filename, &n);
        generate = 0;
        break;
    case 2:
        printf("Enter the size of the matrix (power of 2): ");
        scanf("%u", &n);
        if (!isPowerOfTwo(n))
        {
            printf("Matrix isn't a power of 2, use courageux.exe\n");
            return 0;
        }
        A = dataMatrix(n);
        for (uint32_t i = 0; i < n; i++)
        {
            for (uint32_t j = 0; j < n; j++)
            {
                A[i][j] = ((double)rand() / RAND_MAX) * 20.0 - 10.0;
            }
        }
        generate = 1;
        break;
    default:
        printf("Invalid choice. Exiting.\n");
        return 0;
    }

    double **AT = transpose(A, n);
    double **ATA = multiply(AT, A, n);

    start = clock();
    double **ATA_1 = inverse(ATA, n);
    end = clock();

    double **A_1 = multiply(ATA_1, AT, n);
    double **result = multiply(A_1, A, n);

    printf("1. Print to file\n");
    printf("2. Display in terminal\n");
    printf("Enter your choice (1 or 2): ");

    int outputChoice;
    scanf("%d", &outputChoice);

    switch (outputChoice)
    {
    case 1:
    {
        char inverse_filename[50];
        char identity_filename[50];
        snprintf(inverse_filename, sizeof(inverse_filename), "./results/inverse_%u.txt", n);
        snprintf(identity_filename, sizeof(identity_filename), "./results/identity_%u.txt", n);
        printMatrixToFile(inverse_filename, A_1, n);
        printMatrixToFile(identity_filename, result, n);

        if (generate == 1)
        {
            char generate_filename[50];
            snprintf(generate_filename, sizeof(generate_filename), "./results/generated_%u.txt", n);
            printMatrixToFile(generate_filename, A, n);
        }

        break;
    }
    case 2:
        if (n <= 8)
        {
            printf("Original Matrix :\n");
            printMatrix(A, n);
            printf("Inverted Matrix :\n");
            printMatrix(A_1, n);
            printf("Identity Matrix :\n");
            printMatrix(result, n);
        }
        else
        {
            printf("Original Matrix :\n");
            printMatrixFormatted(A, n);
            printf("Inverted Matrix :\n");
            printMatrixFormatted(A_1, n);
            printf("Identity Matrix :\n");
            printMatrixFormatted(result, n);
        }
        break;
    default:
        printf("Invalid output choice. Exiting.\n");
        break;
    }

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Inversion duration  : %f seconds\n", cpu_time_used);

    freeMatrix(A, n);
    freeMatrix(AT, n);
    freeMatrix(ATA, n);
    freeMatrix(ATA_1, n);
    freeMatrix(A_1, n);
    freeMatrix(result, n);

    return 0;
}