// Ax = b

#include <unistd.h>
#include "lib_poisson1D.h"
#include "blaslapack_headers.h"

#define VERBOSE 1

void print_ge_matrix(const unsigned m, const unsigned n, double *A)
{
  for (size_t i = 0; i < m; i++)
  {
    for (size_t j = 0; j < n; j++)
    {
      printf("%9lf ", A[i * n + j]);
    }
    printf("\n");
  }
}

void print_gb_matrix(const unsigned m, const unsigned lab, double *A)
{
  for (size_t i = 0; i < lab; i++)
  {
    for (size_t j = 0; j < m; j++)
    {
      printf("%9lf ", A[j * lab + i]);
    }
    printf("\n");
  }
}

void tril(const unsigned m, const double *A, double *Alower, const unsigned lab, const unsigned kl, const unsigned ku)
{
  for (unsigned i = 0; i < m * m; i++)
    Alower[i] = 0.0;

  for (unsigned i = 0; i < kl + 1; i++)
    for (unsigned j = 0; j < m - i; j++)
      Alower[(j + i) * m + j] = A[(j * lab) + ku + i];
}

void iter_richardson_gb(const unsigned m, const unsigned lab, const unsigned kl, const unsigned ku, double *A, double *x, double *b, double omega, double *tmp)
{
  // x = x + omega * (b - A * x);

  // tmp = -A*x
  cblas_dgbmv(CblasColMajor, CblasNoTrans, m, m, kl, ku, -1.0, A, lab, x, 1, 0.0, tmp, 1);

  // b + tmp
  cblas_daxpy(m, 1.0, b, 1, tmp, 1);

  // x + omega * tmp
  cblas_daxpy(m, omega, tmp, 1, x, 1);
}

void iter_GaussSeidel_gb(const unsigned m, const unsigned lab, const unsigned kl, const unsigned ku, double *A, double *x, double *b, double *Minv, double *tmp)
{
  // x = x + M⁻¹ * (b - A * x);

  // tmp = -A*x
  cblas_dgbmv(CblasColMajor, CblasNoTrans, m, m, kl, ku, -1.0, A, lab, x, 1, 0.0, tmp, 1);

  // tmp = b + tmp
  cblas_daxpy(m, 1.0, b, 1, tmp, 1);

  double buffer[m];
  // buffer = M⁻¹ * (b - A*x)
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, m, 1.0, Minv, m, tmp, 1, 0.0, buffer, 1);

  // x + buffer
  cblas_daxpy(m, 1.0, buffer, 1, x, 1);
}

const char check_result(const double *bex, const double *b, const unsigned m, const double error)
{
  for (unsigned i = 0; i < m; i++)
  {
    if (fabs(b[i] - bex[i]) > error)
    {
      //printf("Error: %lf %lf\n", b[i], bex[i]);
      return 0;
    }
  }
  return 1;
}

int main(void)
{
  printf("Ax = b\n");
  printf("Richardson : x = x + omega * (b - A * x)\n");
  printf("Jacobi, Gauss Seidel : x = x + M^(-1) * (b - A * x)\n");
  srand(getpid());
  unsigned m = 5;                  // taille de la matrice m * m
  unsigned kl = 1;                 // nombre de sousdiagonal inferieur
  unsigned ku = 1;                 // nombre de sousdiagonal superieur
  unsigned kv = 0;                 // nombre de sousdiagonal en plus
  unsigned lab = kl + ku + kv + 1; // nombre de sousdiagonale

  double *AB_j = (double *)malloc(sizeof(double) * m * lab);
  double *AB_gs = (double *)malloc(sizeof(double) * m * lab);

  set_GB_operator_colMajor_poisson1D(AB_j, &lab, &m, &kv);
  set_GB_operator_colMajor_poisson1D(AB_gs, &lab, &m, &kv);

  if (VERBOSE)
  {
    printf ("\n\033[35;01mInitialisation\033[00m\n");
    printf("Matrix Band Jacobi: \n");
    print_gb_matrix(m, lab, AB_j);
    printf("\nMatrix Band Gauss Seidel: \n");
    print_gb_matrix(m, lab, AB_gs);
  }

  double *b = (double *)malloc(sizeof(double) * m);
  double *x_gb_j = (double *)malloc(sizeof(double) * m);
  double *x_gb_gs = (double *)malloc(sizeof(double) * m);
  double *tmp_gb_j = (double *)malloc(sizeof(double) * m);
  double *tmp_gb_gs = (double *)malloc(sizeof(double) * m);

  double *Minv = (double *)malloc(sizeof(double) * m * m);

  for (unsigned i = 0; i < m; i++)
  {
    b[i] = (double)rand() / (double)RAND_MAX;
    x_gb_j[i] = x_gb_gs[i] = 0.0;
  }

  const double omega = 2.0 / (eigmin_poisson1D(&m) + eigmax_poisson1D(&m));
  tril(m, AB_gs, Minv, lab, kl, ku);
  LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'U', 'N', m, Minv, m);

  if (VERBOSE)
  {
    printf("\nb = \n");
    print_ge_matrix(m, 1, b);
    printf("\nomega = %lf\n", omega);
    printf("\nM^(-1) = \n");
    print_ge_matrix(m, m, Minv);
  }

  FILE *f_j = fopen("jacobi.txt", "w");
  for (int i = 0; i < 300; i++)
  {
    iter_richardson_gb(m, lab, kl, ku, AB_j, x_gb_j, b, 0.5, tmp_gb_j);
    const double norme_tmp_j = cblas_dnrm2(m, tmp_gb_j, 1);
    fprintf(f_j, "%d %lf\n", i, norme_tmp_j);
  }
  fclose(f_j);

  FILE *f_gs = fopen("gaussSeidel.txt", "w");
  for (int i = 0; i < 100; i++)
  {
    iter_GaussSeidel_gb(m, lab, kl, ku, AB_gs, x_gb_gs, b, Minv, tmp_gb_gs);
    const double norme_tmp_gs = cblas_dnrm2(m, tmp_gb_gs, 1);
    fprintf(f_gs, "%d %lf\n", i, norme_tmp_gs);
  }
  fclose(f_gs);

  if (VERBOSE)
  {
    printf ("\n\033[35;01mRésultat\033[00m");
    printf("\nx_gb_j = \n");
    print_ge_matrix(m, 1, x_gb_j);
    printf("\nx_gb_gs = \n");
    print_ge_matrix(m, 1, x_gb_gs);
    printf("\n");
  }

  printf("--------- Test -----------\n");

  if (VERBOSE)
    printf("A * ~x = ~b | b = ~b ??\n");

  double *res_gb_j = (double *)malloc(sizeof(double) * m);
  double *res_gb_gs = (double *)malloc(sizeof(double) * m);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, m, m, kl, ku, 1.0, AB_j, lab, x_gb_j, 1, 0.0, res_gb_j, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, m, m, kl, ku, 1.0, AB_gs, lab, x_gb_gs, 1, 0.0, res_gb_gs, 1);

  const double norme_bex = cblas_dnrm2(m, b, 1);
  const double norme_bgb_j = cblas_dnrm2(m, res_gb_j, 1);
  const double norme_bgb_gs = cblas_dnrm2(m, res_gb_gs, 1);

  printf("Norme b exacte : %f\n\n", norme_bex);

  printf("Jacobi norme(~b): %lf\n", norme_bgb_gs);
  if (check_result(res_gb_j, b, m, 1e-12))
    printf("b et ~b sont egaux\n");
  else
    printf("b et ~b ne sont pas egaux\n");

  printf("\nGauss Seidel norme(~b): %lf\n", norme_bgb_gs);
  if (check_result(res_gb_gs, b, m, 1e-12))
    printf("b et ~b sont egaux\n");
  else
    printf("b et ~b ne sont pas egaux\n");

  if (VERBOSE)
  {
    printf("\n~b = \n");
    print_ge_matrix(m, 1, res_gb_j);
    printf("\n");
    print_ge_matrix(m, 1, res_gb_gs);
  }

  printf("\n--------- End -----------\n");

  free(AB_j);
  free(AB_gs);
  free(b);
  free(x_gb_j);
  free(x_gb_gs);
  free(tmp_gb_j);
  free(tmp_gb_gs);
  free(Minv);
  free(res_gb_j);
  free(res_gb_gs);

  return 0;
}