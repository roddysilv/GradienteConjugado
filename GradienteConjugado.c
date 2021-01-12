/*
Compilar :  gcc -o GradienteConjugado GradienteConjugado.c -lm
Uso:        ./GradienteConjugado n min max precisao
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

void salvar_csv(double **m, double *b, double *x, int n)
{
  FILE *A = fopen("A.csv", "w");
  FILE *B = fopen("b.csv", "w");
  FILE *X = fopen("x.csv", "w");

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (j < n - 1)
      {
        fprintf(A, "%f,", m[i][j]);
      }
      else
      {
        fprintf(A, "%f", m[i][j]);
      }
    }
    fprintf(A, "\n");
  }

  for (int i = 0; i < n; i++)
  {
    if (i < n - 1)
    {
      fprintf(B, "%f,", b[i]);
      fprintf(X, "%f,", x[i]);
    }
    else
    {
      fprintf(B, "%f", b[i]);
      fprintf(X, "%f", x[i]);
    }
  }

  fclose(A);
  fclose(B);
  fclose(X);
}


double produtoVetorial(double *a, double *b, int n)
{
  int i;
  double r = 0;
  for (i = 0; i < n; i++)
  {
    r += a[i] * b[i];
  }
  return r;
}

void transposta(double **a, int n)
{
  int i, j;
  double temp;

  for (i = 0; i < n; i++)
  {
    for (j = i + 1; j < n; j++)
    {
      temp = a[i][j];
      a[i][j] = a[j][i];
      a[j][i] = temp;
    }
  }
}

void produtoMatrizVetor(double **A, double *b, double *c, int n)
{
  int i, j;
  double r;
  for (i = 0; i < n; i++)
  {
    r = 0;
    for (j = 0; j < n; j++)
    {
      r += A[i][j] * b[j];
    }
    c[i] = r;
  }
}

void mspd(double **a, int n)
{
  double **b = malloc(n * sizeof(double *));
  double **id = malloc(n * sizeof(double *));

  for (int i = 0; i < n; i++)
  {
    b[i] = malloc(n * sizeof(double));
    id[i] = malloc(n * sizeof(double));
  }

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      b[i][j] = a[i][j];
      if (i == j)
      {
        id[i][j] = 1;
      }
      else
      {
        id[i][j] = 0;
      }
    }
  }

  transposta(b, n);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      a[i][j] = a[i][j] * b[i][j] - id[i][j];
    }
  }

  for (int i = 0; i < n; i++)
  {
    free(b[i]);
    free(id[i]);
  }

  free(b);
  free(id);
}

void gradienteConjugado(double **A, double *b, double *x, int n, double precisao)
{
  int i, j;

  double *r, *p, *c, *Ap;

  r = (double *)malloc(n * sizeof(double));
  p = (double *)malloc(n * sizeof(double));
  c = (double *)malloc(n * sizeof(double));
  Ap = (double *)malloc(n * sizeof(double));

  produtoMatrizVetor(A, x, c, n);

  for (i = 0; i < n; i++)
  {
    r[i] = b[i] - c[i];
    p[i] = r[i];
  }

  double rsold = produtoVetorial(r, r, n);

  j = 0;

  while (sqrt(produtoVetorial(r, r, n)) > precisao)
  {

    j++;

    produtoMatrizVetor(A, p, Ap, n);

    double alpha = rsold / (produtoVetorial(p, Ap, n));

    for (i = 0; i < n; i++)
    {
      x[i] = x[i] + alpha * p[i];

      r[i] = r[i] - alpha * Ap[i];
    }

    double rsnew = produtoVetorial(r, r, n);

    for (i = 0; i < n; i++)
    {
      p[i] = r[i] + (rsnew / rsold) * p[i];
    }

    rsold = rsnew;
  }

  free(Ap);
  free(r);
  free(p);
  free(c);

  printf("Iterações: %d\n", j);
}

int main(int argc, char *argv[])
{
  clock_t start, end;

  double cpu_time_used;

  int n = strtol(argv[1], NULL, 10);

  double precisao = strtol(argv[4], NULL, 10);

  int i, j;

  //srand(time(0));

  double **m = malloc(n * sizeof(double *));

  for (i = 0; i < n; i++)
  {
    m[i] = malloc(n * sizeof(double));
  }

  double *x = (double *)malloc(n * sizeof(double));
  double *b = (double *)malloc(n * sizeof(double));

  int min, max;

  min = strtol(argv[2], NULL, 10);

  max = strtol(argv[3], NULL, 10);

  for (i = 0; i < n; i++)
  {
    b[i] = min + ((double)rand() / (double)(RAND_MAX)) * (max - (min));
    x[i] = 0;
    for (j = 0; j < n; j++)
    {
      m[i][j] = min + ((double)rand() / (double)(RAND_MAX)) * (max - (min));
    }
  }

  mspd(m, n);

  start = clock();

  gradienteConjugado(m, b, x, n, precisao);

  end = clock();

  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

  printf("t = %f s\n", cpu_time_used);

  salvar_csv(m, b, x, n);

  for (i = 0; i < n; i++)
  {
    free(m[i]);
  }

  free(x);
  free(b);
  free(m);

  return 0;
}