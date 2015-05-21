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
      printf("%.15g", M(a, i, k, n));
      if (k != n - 1)
        printf(" ");
    }
    printf("\n");
  }
}

// Prints a vector
void print_vector(double *a, int n) {
  int i;
  for (i = 0; i < n; i++) {
    printf("%.15g\n", a[i]);
  }
}

// Copy a int b
void copy_vector(double *a, double *b, int n) {
  int i;
  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

// Copy a into b
void copy_matrix(double *a, double *b, int n) {
  int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      M(b, i, k, n) = M(a, i, k, n);
}

// Zeroes the matrix
void zero_matrix(double *a, int n) {
  int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      M(a, i, k, n) = 0;
}

// Transposes a into b
void transpose_matrix(double *a, double *b, int n) {
  double *temp = malloc(sizeof(double) * n * n);
  int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      M(temp, i, k, n) = M(a, k, i, n);
  copy_matrix(temp, b, n);
  free(temp);
}

// Multiplies matrix by scalar
void multiply_scalar(double *a, double s, int n, double *res) {
  int i, k;
  zero_matrix(res, n);
  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++) {
      M(res, i, k, n) = M(a, i, k, n) * s;
    }
  }
}

// Adds a + b = res
void add_matrices(double *a, double *b, int n, double *res) {
  int i, k;
  zero_matrix(res, n);
  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++) {
      M(res, i, k, n) = M(a, i, k, n) + M(b, i, k, n);
    }
  }
}

// Multiplies a * b = res
void multiply_matrix(double *a, double *b, double *res, int n) {
  int i, k, r;
  zero_matrix(res, n);
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      for (r = 0; r < n; r++)
        M(res, i, k, n) += M(a, i, r, n) * M(b, r, k, n);
}

void multiply_matrix_by_vector(double *a, double *b, double *res, int n) {
  int i, k;
  for (i = 0; i < n; i++) {
    res[i] = 0;
  } 
  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++) {
      res[i] += b[k] * M(a, i, k, n);
    }
  }
}

// Makes id into identity
void identity_matrix(double *id, int n) {
  int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      M(id, i, k, n) = (i == k) ? 1 : 0;
}

// Checks if float is equal to number to the absolute error e
char absisequal(double a, double b, double e) {
  return fabs(b - a) < e;
}

// Checks if float is "equal" to number to the relative error re
char relisequal(double a, double b, double re) {
  return fabs((b - a) / b) < re;
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

void backward_euler(double *a, int n, double *y0, double t, int m, double *y) {
  int i;
  double step = t/m;
  double *temp1 = malloc(sizeof(double) * n * n);
  double *temp2 = malloc(sizeof(double) * n * n);
  double *temp3 = malloc(sizeof(double) * n * n);
  double *temp4 = malloc(sizeof(double) * n);
  double *temp5 = malloc(sizeof(double) * n);
  double *id = malloc(sizeof(double) * n * n);

  identity_matrix(id, n);

  multiply_scalar(a, step, n, temp1);
  add_matrices(temp1, id, n, temp2); 
  inv_double_gs(temp2, n, temp3, temp1);
  multiply_matrix_by_vector(temp1, y0, temp4, n);

  for (i = 1; i < m; i++) {
    multiply_matrix_by_vector(temp1, temp4, temp5, n); 
    copy_vector(temp5, temp4, n);
  }

  copy_vector(temp4, y, n);

  free(id);
  free(temp1);
  free(temp2);
  free(temp3);
  free(temp4);
  free(temp5);
}  


void stiff_solve(double *a, int n, double *y0, double t, int m, double *y) {
  int i;
  // Richardson extrapolation.
  double *t1 = malloc(sizeof(double) * n);
  double *t2 = malloc(sizeof(double) * n);
  backward_euler(a, n, y0, t, m, t1);
  backward_euler(a, n, y0, t, m * 2, t2);
  for (i = 0; i < n; i++) {
    t2[i] *= 2;
  }
  for (i = 0; i < n; i++) {
    y[i] = t2[i] - t1[i];
  }
  free(t1);
  free(t2);
}
