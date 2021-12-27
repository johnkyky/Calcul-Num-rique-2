#include "lib_poisson1D.h"
#include "blaslapack_headers.h"


//genere dans AB une matrice tridiagonal de taille (n * n) au format general
//genere dans AB une matrice tridiagonal de taille (n * n) au format bande
void generate_band_matrix(unsigned n, double *AB, double *AB_bande)
{
    for (unsigned i = 0; i < n * n; i++)
    {
        AB_bande[i] = 0.0;
    }

    for (unsigned i = 0; i < n - 1; i++)
    {
        AB_bande[i * n + (i + 1)] = (double)rand() / (double)RAND_MAX;
        AB_bande[i * n + i] = (double)rand() / (double)RAND_MAX;
        AB_bande[(i + 1) * n + i] = (double)rand() / (double)RAND_MAX;
    }
    AB[(n - 1) * 3 + 1] = AB_bande[n * n - 1] = (double)rand() / (double)RAND_MAX;

    AB[0] = 0.0;
    for(unsigned i = 0; i < n - 1; i++)
    {
        AB[(i + 1) * 3] = AB_bande[i * n + (i + 1)];
        AB[i * 3 + 1] = AB_bande[i * n + i];
        AB[i * 3 + 2] = AB_bande[(i + 1) * n + i];
    }

    printf("Matrice general\n");
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            printf("%f ", AB_bande[i * n + j]);
        printf("\n");
    }

    printf("\nMatrice general bande (col major)\n");
    for (unsigned i = 0; i < 3 * n; i++)
        printf("%f\n", AB[i]);
}

int main(int argc, char *argv[])
{
    int nbpoints = 0, la = 0;
    int ku = 0, kl = 0, kv = 0, lab = 0;

    nbpoints = 102;
    la = nbpoints - 97;

    printf("--------- DGBMV ---------\n\n");

    kv = 0;
    ku = 1;
    kl = 1;
    lab = kv + kl + ku + 1;

    double *AB = (double *)malloc(sizeof(double) * lab * la);
    double *x = (double *)malloc(sizeof(double) * la);
    double *y = (double *)malloc(sizeof(double) * la);

    double *AB_ex = (double *)malloc(sizeof(double) * la * la);
    double *x_ex = (double *)malloc(sizeof(double) * la);
    double *y_ex = (double *)malloc(sizeof(double) * la);

    generate_band_matrix(la, AB, AB_ex);

    for (unsigned i = 0; i < la; i++)
    {
        x_ex[i] = x[i] = (double)rand() / (double)RAND_MAX;
        y_ex[i] = y[i] = (double)rand() / (double)RAND_MAX;
    }

    printf("\n");
    for(int i = 0; i < la; i++)
        printf("%f\n", x_ex[i]);
    printf("\n");
    for(int i = 0; i < la; i++)
        printf("%f\n", y_ex[i]);
    
    printf("\n");
    for(int i = 0; i < la; i++)
        printf("%f\n", x[i]);
    printf("\n");
    for(int i = 0; i < la; i++)
        printf("%f\n", y[i]);

    double alpha = (double)rand() / (double)RAND_MAX;
    double beta = (double)rand() / (double)RAND_MAX;

    int row = 0;
    if (row == 1)
    { // LAPACK_ROW_MAJOR
        cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, alpha, AB, lab, x, 1, beta, y, 1);
    }
    else
    { // LAPACK_COL_MAJOR
        cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB, lab, x, 1, beta, y, 1);
    }


    printf("\n\n--------- TEST -----------\n");

    cblas_dgemv(CblasRowMajor, CblasNoTrans, la, la, alpha, AB_ex, la, x_ex, 1, beta, y_ex, 1);

    char is_egal = 1;
    for (unsigned i = 0; i < la; i++)
    {
        if (fabs(y_ex[i] - y[i]) > 1e-15)
            is_egal = 0;
    }

    if(is_egal)
        printf("Les vecteurs y et y_exact sont egaux\n");
    else
        printf("Les vecteurs y et y_exact sont pas egaux\n");

    double norme_y = cblas_dnrm2(la, y, 1);
    double norme_yex = cblas_dnrm2(la, y_ex, 1);

    printf("Norme y : %f\n", norme_y);
    printf("Norme y_ex : %f\n", norme_yex);

    free(AB);
    free(x);
    free(y);
    free(AB_ex);
    free(x_ex);
    free(y_ex);

    printf("\n--------- End -----------\n");

    return 0;
}
