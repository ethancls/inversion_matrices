#include "mat.h"

typedef double **(*multFunc)(double **, double **, uint32_t);

// This function computes the reciprocal (1/x) using magic value **NOT USED HERE**
double reciprocal(double x)
{
    // Convert the bits of x to an integer, and subtract it from a constant magic number.
    uint64_t u = 0x7FDE623822FC16E6 - (*((uint64_t *)&x));

    // Convert the bits back to a double precision floating-point number.
    double f = *((double *)&u);

    // Apply three iterations of the Newton-Raphson method to refine the reciprocal estimate.
    for (int i = 0; i < 3; i++)
    {
        f = 2.0 * f - f * f * x;
    }

    return f;
}

// Pad matrix with identity to the nearest power of 2, bcs it dont change the determinant of the matrix
double **padMatrix(double **A, uint32_t originalN, uint32_t *newN)
{
    uint32_t newSize = 1;
    while (newSize < originalN)
    {
        newSize *= 2;
    }

    // Set the new size
    *newN = newSize;

    // Create a new matrix of the padded size
    double **result = dataMatrix(newSize);

    // Fill the remaining space with the identity matrix
    for (uint32_t i = 0; i < newSize; i++)
    {
        for (uint32_t j = 0; j < newSize; j++)
        {
            result[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Copy the original matrix into the top-left corner
    for (uint32_t i = 0; i < originalN; i++)
    {
        for (uint32_t j = 0; j < originalN; j++)
        {
            result[i][j] = A[i][j];
        }
    }

    return result;
}

// Multiply matrix with Strassen algorithm using 7 multiplications
double **strassen(double **A, double **B, uint32_t size)
{
    // Base case : 1x1 matrix
    if (size == 1)
    {
        double **C = dataMatrix(1);
        C[0][0] = A[0][0] * B[0][0];
        return C;
    }

    // Split matrix
    double **A11, **A12, **A21, **A22;
    double **B11, **B12, **B21, **B22;

    A11 = dataMatrix(size / 2);
    A12 = dataMatrix(size / 2);
    A21 = dataMatrix(size / 2);
    A22 = dataMatrix(size / 2);

    B11 = dataMatrix(size / 2);
    B12 = dataMatrix(size / 2);
    B21 = dataMatrix(size / 2);
    B22 = dataMatrix(size / 2);

    for (uint32_t i = 0; i < size / 2; i++)
    {
        for (uint32_t j = 0; j < size / 2; j++)
        {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + size / 2];
            A21[i][j] = A[i + size / 2][j];
            A22[i][j] = A[i + size / 2][j + size / 2];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + size / 2];
            B21[i][j] = B[i + size / 2][j];
            B22[i][j] = B[i + size / 2][j + size / 2];
        }
    }

    double **t1 = add(A11, A22, size / 2);
    double **t2 = add(B11, B22, size / 2);
    double **M1 = strassen(t1, t2, size / 2);

    double **t3 = add(A21, A22, size / 2);
    double **M2 = strassen(t3, B11, size / 2);

    double **t4 = subtract(B12, B22, size / 2);
    double **M3 = strassen(A11, t4, size / 2);

    double **t5 = subtract(B21, B11, size / 2);
    double **M4 = strassen(A22, t5, size / 2);

    double **t6 = add(A11, A12, size / 2);
    double **M5 = strassen(t6, B22, size / 2);

    double **t7 = subtract(A21, A11, size / 2);
    double **t8 = add(B11, B12, size / 2);
    double **M6 = strassen(t7, t8, size / 2);

    double **t9 = subtract(A12, A22, size / 2);
    double **t10 = add(B21, B22, size / 2);
    double **M7 = strassen(t9, t10, size / 2);

    // Calculate submatrices of the result
    double **t11 = add(M1, M4, size / 2);
    double **t12 = subtract(t11, M5, size / 2);
    double **C11 = add(t12, M7, size / 2);

    double **C12 = add(M3, M5, size / 2);
    double **C21 = add(M2, M4, size / 2);

    double **t13 = subtract(M1, M2, size / 2);
    double **t14 = add(t13, M3, size / 2);
    double **C22 = add(t14, M6, size / 2);

    freeMatrix(A11, size / 2);
    freeMatrix(A12, size / 2);
    freeMatrix(A21, size / 2);
    freeMatrix(A22, size / 2);

    freeMatrix(B11, size / 2);
    freeMatrix(B12, size / 2);
    freeMatrix(B21, size / 2);
    freeMatrix(B22, size / 2);

    freeMatrix(t1, size / 2);
    freeMatrix(t2, size / 2);
    freeMatrix(t3, size / 2);
    freeMatrix(t4, size / 2);
    freeMatrix(t5, size / 2);
    freeMatrix(t6, size / 2);
    freeMatrix(t7, size / 2);
    freeMatrix(t8, size / 2);
    freeMatrix(t9, size / 2);
    freeMatrix(t10, size / 2);
    freeMatrix(t11, size / 2);
    freeMatrix(t12, size / 2);
    freeMatrix(t13, size / 2);
    freeMatrix(t14, size / 2);

    freeMatrix(M1, size / 2);
    freeMatrix(M2, size / 2);
    freeMatrix(M3, size / 2);
    freeMatrix(M4, size / 2);
    freeMatrix(M5, size / 2);
    freeMatrix(M6, size / 2);
    freeMatrix(M7, size / 2);

    double **result = dataMatrix(size);
    combineMatrices(result, C11, C12, C21, C22, size);

    freeMatrix(C11, size / 2);
    freeMatrix(C12, size / 2);
    freeMatrix(C21, size / 2);
    freeMatrix(C22, size / 2);

    return result;
}

// Inverse with pointer to mult function
double **inverse(double **A, uint32_t size, multFunc mult)
{
    // Base case : 1x1 matrix
    if (size == 1)
    {
        if (A[0][0] == 0)
        {
            printf("Non inversible !");
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
        double **B_1 = inverse(B, size / 2, mult);

        // CB_1, B_1CT, CB_1CT, S
        double **CB_1 = mult(C, B_1, size / 2);
        double **B_1CT = transpose(CB_1, size / 2);
        double **CB_1CT = mult(C, B_1CT, size / 2);
        double **S = subtract(D, CB_1CT, size / 2);

        // Recursive inversion of S
        double **S_1 = inverse(S, size / 2, mult);

        // S_1CB_1, B_1CTS_1, B_1CTS_1CB_1
        double **S_1CB_1 = mult(S_1, CB_1, size / 2);
        double **B_1CTS_1 = transpose(S_1CB_1, size / 2);
        double **B_1CTS_1CB_1 = mult(B_1CTS_1, CB_1, size / 2);

        // B−1 + B−1CT S−1CB−1
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

int main()
{
    srand(time(NULL));

    uint32_t n, paddedSize;
    clock_t start, end;
    double cpu_time_used;
    int choice, generate;
    char filename[100];
    double **A;

    printf("\n1. Use existing matrix file\n");
    printf("2. Generate a matrix\n");
    printf("Enter your choice (1 or 2): ");
    scanf("%d", &choice);
    printf("\n--------------------------------\n");

    switch (choice)
    {
    case 1:
        printf("\nEnter the name of the matrix file: ");
        scanf("%s", filename);
        A = readFromFile(filename, &n);
        generate = 0;
        break;
    case 2:
        printf("\nEnter the size of the matrix: ");
        scanf("%u", &n);
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

    if (!isPowerOfTwo(n))
    {
        printf("\nMatrix size is not a power of 2 -- padding with identity matrix\n");
        double **paddedMatrix = padMatrix(A, n, &paddedSize);
        free(A);
        A = paddedMatrix;
    }
    else
    {
        paddedSize = n;
    }

    int funChoice;
    printf("\n1. Use Strassen algorithm\n");
    printf("2. Use standard multiplication\n");
    printf("Enter your choice (1 or 2): ");
    scanf("%d", &funChoice);
    printf("\n--------------------------------\n");

    double **AT, **ATA, **ATA_1, **A_1, **I;

    switch (funChoice)
    {
    case 1:
    {
        multFunc mult = strassen;
        AT = transpose(A, paddedSize);
        ATA = strassen(AT, A, paddedSize);

        start = clock();
        ATA_1 = inverse(ATA, paddedSize, mult);
        end = clock();

        A_1 = strassen(ATA_1, AT, paddedSize);
        I = strassen(A_1, A, paddedSize);
    }
    break;
    case 2:
    {
        multFunc mult = multiply;
        AT = transpose(A, paddedSize);
        ATA = multiply(AT, A, paddedSize);

        start = clock();
        ATA_1 = inverse(ATA, paddedSize, mult);
        end = clock();

        A_1 = multiply(ATA_1, AT, paddedSize);
        I = multiply(A_1, A, paddedSize);
    }
    break;
    default:
        printf("Invalid choice. Exiting.\n");
        break;
    }

    printf("\n1. Print to file\n");
    printf("2. Display in terminal\n");
    printf("Enter your choice (1 or 2): ");

    int outputChoice;
    scanf("%d", &outputChoice);
    printf("\n--------------------------------\n");

    switch (outputChoice)
    {
    case 1:
    {
        char inverse_filename[50];
        char identity_filename[50];
        snprintf(inverse_filename, sizeof(inverse_filename), "./results/inverse_%u.txt", n);
        snprintf(identity_filename, sizeof(identity_filename), "./results/identity_%u.txt", n);
        printMatrixToFile(inverse_filename, A_1, n);
        printMatrixToFile(identity_filename, I, n);

        if (generate == 1)
        {
            char generate_filename[50];
            snprintf(generate_filename, sizeof(generate_filename), "./results/generated_%u.txt", n);
            printMatrixToFile(generate_filename, A, n);
        }
    }
    break;

    case 2:
        if (n <= 8)
        {
            printf("Original Matrix :\n");
            printMatrix(A, n);
            printf("Inverted Matrix :\n");
            printMatrix(A_1, n);
            printf("Identity Matrix :\n");
            printMatrix(I, n);
        }
        else
        {
            printf("Original Matrix :\n");
            printMatrixFormatted(A, n);
            printf("Inverted Matrix :\n");
            printMatrixFormatted(A_1, n);
            printf("Identity Matrix :\n");
            printMatrixFormatted(I, n);
        }
        break;
    default:
        printf("Invalid output choice. Exiting.\n");
        break;
    }

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Inversion duration  : %f seconds\n", cpu_time_used);

    freeMatrix(A, paddedSize);
    freeMatrix(AT, paddedSize);
    freeMatrix(ATA, paddedSize);
    freeMatrix(ATA_1, paddedSize);
    freeMatrix(A_1, paddedSize);
    freeMatrix(I, paddedSize);

    return 0;
}