#include "algorithm.h"
#include <pthread.h>
#include <iostream>
#include <limits>
#include <math.h>



void synchronize(int total_threads)
{
        static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
        static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
        static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
        static int threads_in = 0;
        static int threads_out = 0;

        pthread_mutex_lock(&mutex);

        threads_in++;
        if (threads_in >= total_threads)
        {
                threads_out = 0;
                pthread_cond_broadcast(&condvar_in);
        } else
                while (threads_in < total_threads)
                        pthread_cond_wait(&condvar_in,&mutex);

        threads_out++;
        if (threads_out >= total_threads)
        {
                threads_in = 0;
                pthread_cond_broadcast(&condvar_out);
        } else
                while (threads_out < total_threads)
                        pthread_cond_wait(&condvar_out,&mutex);

        pthread_mutex_unlock(&mutex);
}

int parallel_Rotation_Method(   int len, double *copy_matrix, double *inv_matrix, 
                                int number_thread, int _amountThreads,
                                double *sin_rot_angles_str, double *cos_rot_angles_str){

        int j, k;
        int j1,  i1;
        double cos, sin;
        double storage_1, storage_2, sum;
        double inv_storage_1, inv_storage_2;
        int lines_num;
        int colums_num_start, colums_num_end;
        int colums_num_start_inv, colums_num_end_inv;
        int columns_num_end_op, columns_num_end_op_2;
        double sub_const;
        //int flag_exit = 0;
        int i;
       
        //printf(" = %d", _amountThreads );
        for(j = 0; j < len; j++){

                if(number_thread == 0){

                        for(k = j + 1; k < len; k++){
                                
                                sum = copy_matrix[j*len+j]*copy_matrix[j*len+j] + copy_matrix[k*len+j]*copy_matrix[k*len+j];

                                if(fabs(sum) < std::numeric_limits<double>::epsilon()){
                                    //flag[0] = -1;
                                    return -1;
                                }

                                cos = copy_matrix[j*len+j]/sqrt(sum);
                                sin = -copy_matrix[k*len+j]/sqrt(sum);

                                cos_rot_angles_str[k-1] = cos;
                                sin_rot_angles_str[k-1] = sin;

                         
                                storage_1=copy_matrix[j*len+j]; 
                                storage_2=copy_matrix[k*len+j];  
                                copy_matrix[j*len+j] = cos*storage_1 - sin*storage_2;
                                copy_matrix[k*len+j] = sin*storage_1 + cos*storage_2;
                                
                                /*
                                inv_storage_1=inv_matrix[j*len+j]; 
                                inv_storage_2=inv_matrix[k*len+j];  
                                inv_matrix[j*len+j] = cos*inv_storage_1 - sin*inv_storage_2;
                                inv_matrix[k*len+j] = sin*inv_storage_1 + cos*inv_storage_2;
                                */
                                

                        }

                        //flag[1] = 0;

                }



                synchronize(_amountThreads);

                
                
                /*parallel matrix transformation by columns*/
                
                if(len - j - 1 > _amountThreads){
                        lines_num = j;
                        
                        colums_num_start = 1 + j + ( (len - j)/(_amountThreads) )*number_thread;
                        
                        
                        columns_num_end_op = len;
                        columns_num_end_op_2 = 1 + j + ( (len - j)/(_amountThreads) )*(number_thread + 1);
                        if( columns_num_end_op_2 <= len){
                                columns_num_end_op = columns_num_end_op_2;
                        }
                        if(number_thread == _amountThreads - 1){
                                columns_num_end_op = len;
                        }

                        colums_num_end = columns_num_end_op;
                }


                if(len - j - 1 <= _amountThreads){

                                lines_num = j;
                         
                                if(number_thread <= len - j - 1){
                                        colums_num_start = j + 1 + number_thread;
                                        colums_num_end = j + 2 + number_thread;
                                        if(number_thread == len - j  -1){
                                            colums_num_end = len;
                                        }
                                }
                                if(number_thread > len - j  - 1){
                                        colums_num_start = 1;
                                        colums_num_end = 0;
                                }
                }



                for(j1 = colums_num_start; j1 < colums_num_end; j1++){

                        for(i1 = lines_num + 1; i1 < len; i1++){

                                
                                cos = cos_rot_angles_str[i1-1];
                                sin = sin_rot_angles_str[i1-1];

                                storage_1 = copy_matrix[lines_num*len +j1];
                                storage_2 = copy_matrix[i1*len + j1];

                                        
                                copy_matrix[lines_num*len+j1] = cos*storage_1 - sin*storage_2;
                                copy_matrix[i1*len+j1] = sin*storage_1 + cos*storage_2;

                        }
                }



                /*parallel transformation of inverse matrix by columns*/
                
                
                if(len > _amountThreads){

                        lines_num = j;

                        colums_num_start_inv = (len/_amountThreads)*number_thread;
                        colums_num_end_inv = (len/_amountThreads)*(number_thread + 1);
                        
                        if(number_thread == _amountThreads - 1){
                                colums_num_end_inv = len;
                        }
                }

                if(len <= _amountThreads){

                        lines_num = j;
                        
                        if(number_thread < len){
                                        colums_num_start_inv = number_thread;
                                        colums_num_end_inv =  1 + number_thread;
                                        if(number_thread == len -1){
                                            colums_num_end_inv = len;
                                        }
                                }
                        if(number_thread >= len){
                                colums_num_start_inv = 1;
                                colums_num_end_inv = 0;
                        }

                }

                //printf("%d  %d %d %d  %d\n", len - j, _amountThreads, colums_num_start_inv, colums_num_end_inv, lines_num);


                for(j1 = colums_num_start_inv; j1 < colums_num_end_inv; j1++ ){

                            for(i1 = lines_num + 1; i1< len ; i1++){

                                cos = cos_rot_angles_str[i1-1];
                                sin = sin_rot_angles_str[i1-1];


                                inv_storage_1=inv_matrix[lines_num*len+j1];
                                inv_storage_2=inv_matrix[i1*len+j1];

                                inv_matrix[lines_num*len+j1] = cos*inv_storage_1 - sin*inv_storage_2;
                                inv_matrix[i1*len+j1] = sin*inv_storage_1 + cos*inv_storage_2;

                            }

                }
                
                

                synchronize(_amountThreads);

        }


        for(k = len - 1 ; k >= 0; k--){

            for(j = k - 1; j >= 0 ; j--){
                sub_const = -1*(copy_matrix[j*len + k]/copy_matrix[k*len + k]);

                for(i = number_thread; i < len; i = i + _amountThreads){
                    inv_matrix[j*len+i] = inv_matrix[j*len+i] + inv_matrix[k*len+i]*sub_const;

                }

            }
        }




        return 0;

        }







































                /*
                 *
                 *
                 *
                 *
                 *
                 *
                 *
                 *for(i = len - 1 ; i>=0; i--){
            for(j = k + 1; j < len; j++){

                for(k = number_thread; k < len; k+=_amountThreads){

                    inv_matrix[i*len+k] -= inv_matrix[j*len+k]*copy_matrix[i*len + j];
                }
            }

            synchronize(_amountThreads);
            if(number_thread == 0){
                for(j = 0; j < len; j++){
                    inv_matrix[i*len+j] = inv_matrix[i*len+j]/copy_matrix[i*len + i];
                }
            }
            synchronize(_amountThreads);

        }
                 *
                 *
                 *
                 *
                if(len - j == _amountThreads){    
                        printf("lox \n");
                
                        if(flag[0] == 1){
                                return 2;
                        }
  
                
                for(j1 = j; j1 < len; j1++){
                        
                        for(k1 = j1 + 1; k1 < len; k1++){
                                
                                sum = copy_matrix[j*len+j]*copy_matrix[j*len+j] + copy_matrix[k*len+j]*copy_matrix[k*len+j];
                                        
                                if(fabs(sum) < 100*std::numeric_limits<double>::epsilon()){
                                        return -1;
                                }
        
                                cos = copy_matrix[j1*len+j1]/sqrt(sum);
                                sin = -copy_matrix[k1*len+j1]/sqrt(sum);
        
                                for(i1 = 0; i1 < len; i1++){
        
                                        storage_1=copy_matrix[j1*len+i1]; 
                                        storage_2=copy_matrix[k1*len+i1];
                                        copy_matrix[j1*len+i1] = cos*storage_1 - sin*storage_2;
                                        copy_matrix[k1*len+i1] = sin*storage_1 + cos*storage_2;


                                        inv_storage_1=inv_matrix[j1*len+i1]; 
                                        inv_storage_2=inv_matrix[k1*len+i1];
                                        inv_matrix[j1*len+i1] = cos*inv_storage_1 - sin*inv_storage_2;
                                        inv_matrix[k1*len+i1] = sin*inv_storage_1 + cos*inv_storage_2;

                                }
        
                        }

                }

                flag[0] = 1;
                return 0;
                }
                */

