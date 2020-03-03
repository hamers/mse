
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
 
using namespace std;
 
extern "C" double function_sum_(long *fsize, double* fvec);
extern "C" void subroutine_sum_(long *fsize, double* fvec, double *fsum);
 

int main(int argc, char ** argv)
{
    long i,size;
    double sum;
    double *vec;
 
    
    
    size = 5000;
 
    vec = new double[size];
 
    for(i=0;i<size;i++)
    {
        vec[i]=1.1234532;
    }
 
    // Call a Fortran function
    sum = function_sum_(&size,vec);
 
    cout << "Calling a Fortran function" << endl;
    cout << "============================" << endl;
    cout << "size = " << size << endl;
    cout << "sum = " << sum << endl << endl;
 
 
    // Call a Fortran subroutine
    subroutine_sum_(&size,vec,&sum);
 
    cout << "Calling a Fortran subroutine" << endl;
    cout << "===============================" << endl;
    cout << "size = " << size << endl;
    cout << "sum = " << sum << endl << endl;
 
    delete[] vec;
 
}
