#ifndef _MAIN_ALGORITHM__15
#define _MAIN_ALGORITHM__15

#include "cmath"

/*
n - размерность матрицы
k - номер потока
p - число потоков
A - матрицы
*/

int parallel_Rotation_Method(int len, double *copy_matrix, double *inv_matrix, 
                            int number_thread, int _amountThreads, 
                            double *sin_rot_angles_str, double *cos_rot_angles_str
                            );

void synchronize(int total_threads);


#endif
