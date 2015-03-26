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

// Multiplies a * b = res
void multiply_matrix(double *a, double *b, double *res, int n) {
  int i, k, r;
  zero_matrix(res, n);
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      for (r = 0; r < n; r++)
        M(res, i, k, n) += M(a, i, r, n) * M(b, r, k, n);
}

// Makes id into identity
void identity_matrix(double *id, int n) {
  int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < n; k++)
      M(id, i, k, n) = (i == k) ? 1 : 0;
}

// Quicksort helpers
void swap(double *tosort, int *columns, int a, int b) {
  double i = tosort[a];
  int j = columns[a];
  tosort[a] = tosort[b]; columns[a] = columns[b];
  tosort[b] = i; columns[b] = j;
}

int partition(double *tosort, int *columns, int l, int h) {
  int pivot = l + (h - l) / 2;
  double v = tosort[pivot];
  swap(tosort, columns, pivot, h);
  int toPut = l;
  
  int i;
  for (i = l; i < h; i++)
    if (tosort[i] > v)
      swap(tosort, columns, i, toPut++);

  swap(tosort, columns, toPut, h);
  return toPut;
}

void quicksort(double *tosort, int *columns, int l, int h) {
  if (l < h) {
    int p = partition(tosort, columns, l, h);
    quicksort(tosort, columns, l, p - 1);
    quicksort(tosort, columns, p + 1, h);
  }
}

// Checks if float is equal to number to the absolute error e
char absisequal(double a, double b, double e) {
  return fabs(b - a) < e;
}

// Checks if float is "equal" to number to the relative error re
char relisequal(double a, double b, double re) {
  return fabs((b - a) / b) < re;
}

// Check if all the off diagonals are equal to 0 after each sweep
char issweepcomplete(double *a, int n) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i != j && !absisequal(0.0, M(a, i, j, n), 1e-16))
        return 0;
    }
  }
  return 1;
}

// Checks if matricies are float equals
char matrix_equals(double *a, double *b, int n) {
  int i, k;
  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++) {
      if (!relisequal(M(a, i, k, n) + 1, M(b, i, k, n) + 1, .00001)) {
        return 0;
      }
    }
  }
  return 1;
}

void jacobi(double *a, int n, double *s, double *u, double *v) {
  int p, q, i, k;
  int first = 1;

  identity_matrix(v, n);
  identity_matrix(u, n);

  while (first == 1 || !issweepcomplete(a, n)) {
    for (p = 0; p < n - 1; p++) {
      for (q = p + 1; q < n; q++) {
        // Implicit 2x2 matrix
        double app = M(a, p, p, n);
        double apq = M(a, p, q, n);
        double aqp = M(a, q, p, n);
        double aqq = M(a, q, q, n);

        // Calculate rotation angles
        double at1 = atan((aqp - apq)/(app + aqq));
        double at2 = atan((apq + aqp)/(app - aqq));
        double theta = 0.5 * (at1 + at2);
        double phi = 0.5 * (at2 - at1);

        // Update rows p and q of A and U
        for (i = 0; i < n; i++) {
          double temp = M(a, p, i, n);
          M(a, p, i, n) = M(a, p, i, n) * cos(theta) + M(a, q, i, n) * sin(theta);
          M(a, q, i, n) = M(a, q, i, n) * cos(theta) - temp * sin(theta);
          temp = M(u, p, i, n);
          M(u, p, i, n) = M(u, p, i, n) * cos(theta) + M(u, q, i, n) * sin(theta);
          M(u, q, i, n) = M(u, q, i, n) * cos(theta) - temp * sin(theta);
        }

        // Update columns p and q of A and V
        for (i = 0; i < n; i++) {
          double temp = M(a, i, p, n);
          M(a, i, p, n) = M(a, i, p, n) * cos(phi) + M(a, i, q, n) * sin(phi);
          M(a, i, q, n) = M(a, i, q, n) * cos(phi) - temp * sin(phi);
          temp = M(v, i, p, n);
          M(v, i, p, n) = M(v, i, p, n) * cos(phi) + M(v, i, q, n) * sin(phi);
          M(v, i, q, n) = M(v, i, q, n) * cos(phi) - temp * sin(phi);
        }
      }
    } 
    first++;
  }

  // Transpose U to make it correct
  transpose_matrix(u, u, n);

  // Make S, make all singular values positive.
  for (i = 0; i < n; i++) {
    if (M(a, i, i, n) < 0) {
      M(a, i, i, n) *= -1;
      // Change appropriate columns
      for (k = 0; k < n; k++) {
        M(u, k, i, n) *= -1;
      }
    }
    s[i] = M(a, i, i, n);
  }

  int *columns = malloc(n * sizeof(int));
  for (i = 0; i < n; i++) columns[i] = i;
  quicksort(s, columns, 0, n-1);
  
  // Tmp matrices to copy into
  double *t1 = malloc(n * n * sizeof(double));
  double *t2 = malloc(n * n * sizeof(double));

  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++) {
      M(t1, k, i, n) = M(u, k, columns[i], n);
      M(t2, k, i, n) = M(v, k, columns[i], n);
    }
  }

  copy_matrix(t1, u, n);
  copy_matrix(t2, v, n);

  free(columns);
  free(t1);
  free(t2);
}

/*
int main(int argc, char **argv) {
  int i, k;
  int n = 4;
  double *a = malloc(sizeof(double) * n * n);
  double *s = malloc(sizeof(double) * n);
  double *u = malloc(sizeof(double) * n * n);
  double *v = malloc(sizeof(double) * n * n);
  double *t = malloc(sizeof(double) * n * n);
  double *t2 = malloc(sizeof(double) * n * n);
  double *id = malloc(sizeof(double) * n * n);
  identity_matrix(id, n);  

  // Generate a random matrix
  for (i = 0; i < n * n; i++) {
    double entry = ((double)rand()/(double)RAND_MAX);
    entry -= RAND_MAX / 2;
    a[i] = entry; 
  }

  copy_matrix(a, t2, n);
  jacobi(a, n, s, u, v);
  copy_matrix(t2, a, n);

  // Tests
  // S should be decreasing
  for (i = 0; i < n - 1; i++) {
    if (s[i] < s[i + 1])
      fprintf(stderr, "S should be decreasing.\n");
  }

  // U and V should be orthoganal.
  transpose_matrix(u, t, n);
  multiply_matrix(u, t, t2, n);
  if (!matrix_equals(t2, id, n)) {
    fprintf(stderr, "U should be orthoganal.\n");
    print_matrix(t2, n); 
  }
  transpose_matrix(v, t, n);
  multiply_matrix(v, t, t2, n);
  if (!matrix_equals(t2, id, n)) {
    fprintf(stderr, "V should be orthoganal.\n");
  }

  for (i = 0; i < n; i++)
    M(t2, i, i, n) = s[i];
  multiply_matrix(u, t2, v, n);
  multiply_matrix(v, t, t2, n);

  if (!matrix_equals(t2, a, n)) {
    fprintf(stderr, "USVt should equal a\n");
  }
  free(a); free(s); free(u); free(v); free(t); free(t2); free(id);
  return 0;
}
*/
