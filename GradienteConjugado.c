#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

void imprimeMatriz(float **m, int n)
{
  printf("Matrix: \n");
  int i, j;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      printf("%0.2f ", m[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void imprimeVetor(float *a, int n)
{
  printf("Vetor: \n");
  int i, j;
  for (i = 0; i < n; i++)
  {
    printf("%0.2f ", a[i]);
  }
  printf("\n\n");
}

float produtoVetorial(float *a, float *b, int n)
{
  int i;
  float r = 0;
  for (i = 0; i < n; i++)
  {
    r += a[i] * b[i];
  }
  return r;
}

void transposta(float **a, int n)
{
  int i, j;
  float temp;

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

void produtoMatrizVetor(float **A, float *b, float *c, int n)
{
  int i, j;
  float r;
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

void mspd(float **a, int n)
{
  float **b = malloc(n * sizeof(float *));
  float **id = malloc(n * sizeof(float *));

  for (int i = 0; i < n; i++)
  {
    b[i] = malloc(n * sizeof(float));
    id[i] = malloc(n * sizeof(float));
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
}

void gradienteConjugado(float **A, float *b, float *x, int n)
{
  int i;

  float *r, *p, *c, *Ap;

  r = (float *)malloc(n * sizeof(float));
  p = (float *)malloc(n * sizeof(float));
  c = (float *)malloc(n * sizeof(float));
  Ap = (float *)malloc(n * sizeof(float));

  produtoMatrizVetor(A, x, c, n);

  for (i = 0; i < n; i++)
  {
    r[i] = b[i] - c[i];
    p[i] = r[i];
  }

  float rsold = produtoVetorial(r, r, n);

  while (sqrt(produtoVetorial(r, r, n)) > 0.0001)
  {

    produtoMatrizVetor(A, p, Ap, n);

    float alpha = rsold / (produtoVetorial(p, Ap, n));

    for (i = 0; i < n; i++)
    {
      x[i] = x[i] + alpha * p[i];

      r[i] = r[i] - alpha * Ap[i];
    }

    float rsnew = produtoVetorial(r, r, n);

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
}

int main(void)
{

  int n = 45;

  int i, j;

  srand(time(0));

  float **m = malloc(n * sizeof(float *));
  
  for (i = 0; i < n; i++)
  {
    m[i] = malloc(n * sizeof(float));
  }

  float *x = (float *)malloc(n * sizeof(float));
  float *b = (float *)malloc(n * sizeof(float));

  int min, max;

  min = 0;

  max = 10;

  for (i = 0; i < n; i++)
  {
    b[i] = min + ((float)rand() / (float)(RAND_MAX)) * (max - (min));
    x[i] = 0;
    for (j = 0; j < n; j++)
    {
      m[i][j] = min + ((float)rand() / (float)(RAND_MAX)) * (max - (min));
    }
  }

  mspd(m, n);

  gradienteConjugado(m, b, x, n);

  imprimeMatriz(m, n);

  imprimeVetor(b, n);

  imprimeVetor(x, n);

  for (i = 0; i < n; i++)
  {
    free(m[i]);
  }

  free(x);
  free(b);
  free(m);

  return 0;
}