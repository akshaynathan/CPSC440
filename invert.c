#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Macro to access element in matrix
#define M(a, r, c, n) a[r * n + c]

/// Helper Functions
/// All functions assume square matrix.

// Prints the matrix
void print_matrix(double *a, int n) {
  int i, k;
  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++) {
      printf("%.10g", M(a, i, k, n));
      if (k != n - 1)
        printf(" ");
    }
    printf("\n");
  }
}

// Gives the length for a column vector
double vector_length(double *a, int n, int column) {
  double sum = 0;
  int i;
  for (i = 0; i < n; i++)
    sum += pow(M(a, i, column, n), 2);
  return sqrt(sum);
}

// Gives dot product of two column vectors
double dot_product(double *a, int acol, double *b, int bcol, int n) {
  double ret = 0;
  int i;
  for (i = 0; i < n; i++)
    ret += M(a, i, acol, n) * M(b, i, bcol, n);
  return ret;
}

// Copies column from matrix b to a
void copy_column(double *a, int acol, double *b, int bcol, int n) {
  int i;
  for (i = 0; i < n; i++)
    M(a, i, acol, n) = M(b, i, bcol, n);
}

// Transposes a into b
void transpose_matrix(double *a, double *b, int n) {
  int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      M(b, i, k, n) = M(a, k, i, n);
}

// Zeroes the matrix
void zero_matrix(double *a, int n) {
  int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      M(a, i, k, n) = 0;
}

// Multiplies a * b = res
void multiply_matrix(double *a, double *b, double *res, int n) {
  int i, k, r;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      for (r = 0; r < n; r++)
        M(res, i, k, n) += M(a, i, r, n) * M(b, r, k, n);
}

// Use double gram schmidt to invert a with orthonormal basis u and inverse b
void inv_double_gs(double *a, int n, double *u, double *b) {
  // Zero u and b for safety
  zero_matrix(u, n);
  zero_matrix(b, n);

  // Create identity matrix
  double *identity = calloc(n * n, sizeof(double));
  int i, j, r, l;
  for (i = 0; i < n; i++)
    M(identity, i, i, n) = 1;
  
  // Make columns orthogonal to other columns
  for (i = 0; i < n; i++) {
    copy_column(u, i, a, i, n);
    // Iterate through previous columns
    for (j = i - 1; j >= 0; j--) {
      // "Double" Gram-Schmidt
      for (l = 0; l < 2; l++) {
        // Subtract the projection of this column
        double p = dot_product(u, i, u, j, n);
        p /= dot_product(u, j, u, j, n);
        for (r = 0; r < n; r++) {
          M(u, r, i, n) -= p * M(u, r, j, n);
          // Apply same operations to identity
          M(identity, r, i, n) -= p * M(identity, r, j, n);
        }          
      }
    }        
  }

  // Normalize all columns
  for (i = 0; i < n; i++) {
    double len = vector_length(u, n, i);
    for (r = 0; r < n; r++) {
      M(u, r, i, n) /= len;
      // Apply same operation to identity matrix
      M(identity, r, i, n) /= len;
    }
  }
  
  // Do A' = G x U*
  double *orthTranspose = calloc(n * n, sizeof(double));
  transpose_matrix(u, orthTranspose, n);
  multiply_matrix(identity, orthTranspose, b, n);

  
  // Cleanup memory
  free(identity);
  free(orthTranspose);
}
