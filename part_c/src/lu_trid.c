#include "lib_poisson1D.h"
#include "blaslapack_headers.h"

void genere_matrix(const unsigned n, double *A)
{
  for (unsigned i = 0; i < n; i++)
  {
    A[i] = 1.0;
    A[n + i] = 2.0;
    A[2 * n + i] = 3.0;
  }

  A[0] = 0.0;
  A[n * 3 - 1] = 0.0;

  for (unsigned i = 0; i < n * 3; i++)
  {
    printf("%f\n", A[i]);
  }
}


void lu_trid(unsigned n, double *A)
{
  for (unsigned i = 1; i < n; i++)
  {
    double buffer1 = A[2 * n + i - 1];
    double buffer2 = A[n + i - 1];
    double buffer3 = A[n + i];

    A[2 * n + i - 1] = A[2 * n + i - 1] / A[n + i - 1];
    buffer1 = A[2 * n + i - 1];
    A[n + i] = A[n + i] - A[2 * n + i - 1] * A[i];

    printf("l = %lf / %lf\n", buffer1, buffer2);
    printf("l = %lf - %lf * %lf\n", buffer2, buffer1, A[i]);

    printf("l = %lf  v = %lf\n", A[2 * n + i - 1], A[n + i]);

    printf("\n");
    for (unsigned i = 0; i < n * 3; i++)
    {
      printf("%f\n", A[i]);
    }
  }
}

int main(void)
{
  const unsigned n = 5;
  double *A = (double *)malloc(sizeof(double) * n * 3);

  genere_matrix(n, A);
  lu_trid(n, A);

  printf("\n");
  for (unsigned i = 0; i < n * 3; i++)
  {
    printf("%f\n", A[i]);
  }

  free(A);

  return 0;
}