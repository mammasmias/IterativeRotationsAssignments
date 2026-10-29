/* dgemm_ is a real BLAS symbol. libira defines ira_dgemm_ instead, so
   this file links only when a BLAS is named on the link line. */
#include <stddef.h>
#include <stdio.h>

extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                   double *alpha, double *a, int *lda, double *b, int *ldb,
                   double *beta, double *c, int *ldc, size_t transa_len,
                   size_t transb_len);

void libira_get_version(char *string, long *date);

int main(void) {
  char trans = 'N';
  int one = 1;
  double alpha = 1.0;
  double beta = 0.0;
  double a = 2.0;
  double b = 3.0;
  double c = 0.0;
  char version[8];
  long date = 0;
  int i;

  for (i = 0; i < 8; i++) {
    version[i] = 0;
  }
  dgemm_(&trans, &trans, &one, &one, &one, &alpha, &a, &one, &b, &one, &beta,
         &c, &one, (size_t)1, (size_t)1);
  libira_get_version(version, &date);
  if (c != 6.0 || date == 0) {
    fprintf(stderr, "c=%g date=%ld\n", c, date);
    return 1;
  }
  return 0;
}
