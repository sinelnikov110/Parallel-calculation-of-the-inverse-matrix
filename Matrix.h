#ifndef MATRIX_H     
#define MATRIX_H

#include <iostream>
#include "algorithm.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <pthread.h>


#include <sys/time.h>
#include <sys/resource.h>


using namespace std;

class Matrix;

long int get_time();
long int get_full_time();

typedef struct{
    
    int i;
    long int *_t_full;
    int number_thread;
    int len;
    int _amountThreads;
    int lines_num;
    int colums_num_start;
    int colums_num_end;
    int inv_colums_num_start;
    int inv_colums_num_end;
    int number_task;
    
    int *flag;
    double *cos_rot_angles_str;
    double *sin_rot_angles_str;

    double *copy_matrix;
    double *inv_matrix;
    
} optionalStruct;


class Matrix{
    
    private:
        int _amountThreads;
        int len;
        int *flag;
        double *sin_rot_angles;
        double *cos_rot_angles;
        
        
        double *_matrix;
        pthread_t *_threads;  // массив потоков
        optionalStruct *_args; // массив структур
        long int _t_full; /* астрономическое время работы всего процесса */
        
        
    public:
        
        Matrix();
        Matrix(const Matrix &A);
        Matrix(int n, int __amountThreads);
        Matrix(double *list,int n, int __amountThreads);
        ~Matrix();
        

        friend double Restnormen(const Matrix &A, const Matrix &B);
        int paralel_inverse_matrix(double *inv_matrix, double *copy_matrix);
        int Rotation_Method(double *inv_matrix, double *copy_matrix, int j);
        

        int InMatrix_formulas(int k);
        int InMatrix_file(string nameOfFile);
        friend void *process_function(void *pa);
        friend void *inv_process_function(void *pa);

        long int get_elapsed();
        


        void Print(int k);
};



#endif
