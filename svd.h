 

#ifndef SVD_H
#define SVD_H
 
void svdcmp(double **a, int m, int n, double w[], double **v);

double pythag(const double a, const double b);

void invcmp(double **a, int n, double **v);
void invcmp2(double **a, int n, double **v);
 
#endif