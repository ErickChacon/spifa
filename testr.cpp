// this is test.cpp
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
using namespace std;
// these declarations are needed as I don't think there is a lapack.h
extern "C" int dpotrf_(char
*
uplo, int
*
n, double
*
a, int
*
lda, int
*
info);
extern "C" int dsyrk_(char
*
uplo, char
*
trans, int
*
n, int
*
k,
double
*
alpha, double
*
a, int
*
lda, double
*
beta, double
*
c, int
*
ldc);
// this is a one-line comment
/*
This is
a multi-line
comment.
*/
// compilation:
// g++ -o test test.cpp -I/usr/share/R/include -llapack -lblas
//    -lRmath -lR -O3 -Wall
int main(){
int size = 8000;
int info = 0;
char uplo = 'U';
char trans = 'N';
double alpha = 1.0;
double beta = 0.0;
double
*
x = new double[size
*
size];
double
*
C = new double[size
*
size];
for(int i = 0; i < size
*
size; i++){
x[i] = rnorm(0.0, 1.0);
C[i] = 0.0;
}
cout << "done with rnorm" << endl;
dsyrk_(&uplo, &trans, &size, &size, &alpha, x, &size, &beta, C, &size);
cout << "done crossprod" << endl;
dpotrf_(&uplo,&size,C,&size,&info);
cout << "done chol" << endl;
return 0;
}
