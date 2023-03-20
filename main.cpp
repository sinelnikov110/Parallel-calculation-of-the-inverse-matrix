//Финальная версия 19,12,2022
#include "Matrix.h"
#include "time.h"

using namespace std;





int main(int argc, char * argv[]){

    printf("Hallo!\n");

    

    if(argc != 5 && argc != 6){
        cout << "initial parameter error";
        return -1;
    }


    int n, m ,k, p;
    double residual;
    long int elapsed;

    string file_name;

        n = atoi(argv[1]);
        p = atoi(argv[2]);
        m = atoi(argv[3]);
        k = atoi(argv[4]);

    if(m > n){
        cout << "initial parameter error";
        return -1;
    }


    printf("p = %d \n",p);

    Matrix obj1(n, p);

    double *inv_matrix = new double[n*n];
    double *copy_matrix = new double[n*n];

    if(k == 0){
        if(argc != 6){
            cout << "initial parameter error";
            delete [] inv_matrix;
            delete [] copy_matrix;
            return -1;
        }
        file_name = argv[5];
        if(obj1.InMatrix_file(file_name) == -1){
            cout << "ERROR MATRIX";
            delete [] inv_matrix;
            delete [] copy_matrix;
            return -1;
        }
    } else{

        if(obj1.InMatrix_formulas(k) == -1 || argc != 5){
        cout << "Formula is incorrect";
        delete [] inv_matrix;
        delete [] copy_matrix;
        return -1;
        }
    }



    obj1.Print(m);

    
    //clock_t time;
    //time = clock();

    obj1.paralel_inverse_matrix(inv_matrix, copy_matrix);

    Matrix obj2(inv_matrix, n, p);

    cout << "invers matrix : "<< endl;
    obj2.Print(m);
    cout << endl;
    
    residual =  Restnormen(obj1, obj2);
    elapsed = obj1.get_elapsed();

    printf("%s : residual = %e elapsed = %.2ld s = %d n = %d m = %d p = %d\n", argv[0], residual, elapsed, n, m ,k, p);

    delete [] inv_matrix;
    delete [] copy_matrix;
    return 1;

}


