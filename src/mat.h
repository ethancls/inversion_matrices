#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define isPowerOfTwo(n) ((n > 0) && ((n & (n - 1)) == 0))

// Allocate a 2D array for double matrix
double **dataMatrix(uint32_t n);

// Free a matrix (2D array)
void freeMatrix(double **mat, uint32_t n);

// Print matrix on terminal
void printMatrix(double **A, uint32_t n);

// Print matrix to a file
void printMatrixToFile(const char *f, double **A, uint32_t n);

// Function to print the matrix in the specified format
void printMatrixFormatted(double **A, uint32_t n);

// Add matrix A to B
double **add(double **A, double **B, uint32_t n);

// Subtract matrix B from A
double **subtract(double **A, double **B, uint32_t n);

// Return the negative matrix (same as -1 scalar multiplication)
double **negate(double **A, uint32_t n);

// Transpose matrix A to tA
double **transpose(double **A, uint32_t n);

// Standard multiplication
double **multiply(double **A, double **B, uint32_t n);

// Read a matrixfrom a file
double **readFromFile(const char *fileName, uint32_t *n);

// Combine quadrants into one matrix of 'size'
void combineMatrices(double **result, double **C11, double **C12, double **C21, double **C22, uint32_t size);