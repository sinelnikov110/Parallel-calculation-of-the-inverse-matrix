#include "Matrix.h"
#include "algorithm.h"
#include <limits>

/*_____________________________________*/
Matrix::Matrix(){
    int i;
    _matrix = new double[4];
    _matrix[0] = 1.;
    _matrix[1] = 0.;
    _matrix[2] = 0.;
    _matrix[3] = 1.;
    len = 4;

    _amountThreads = 1;
    _threads = new pthread_t[_amountThreads];
    _args = new optionalStruct[_amountThreads];
    cos_rot_angles = new double[len - 1];
    sin_rot_angles = new double[len - 1];
    //flag = new int [2];
    //flag[0] = 0;
    //flag[1] = 0;

    for(i = 0; i < _amountThreads; i++){
        _args[i].cos_rot_angles_str = cos_rot_angles;
        _args[i].sin_rot_angles_str = sin_rot_angles;
        //_args[i].flag = flag;
        _args[i].number_task = i;
        _args[i].number_thread = i;
        _args[i]._t_full = &_t_full;
    }



}

/*_____________________________________*/
Matrix::Matrix(const Matrix &A){

        int i, j;

        _amountThreads = A._amountThreads;
        len = A.len;
        
        _matrix = new double[len*len];
        for(i = 0; i < len; i++){
            for(j = 0; j < len; j++){
                _matrix[i*len + j] = A._matrix[i*len+j];
            }
        }
        _threads = new pthread_t[_amountThreads];  // массив потоков
        _args = new optionalStruct[_amountThreads]; // массив структур

        cos_rot_angles = new double[len - 1];
        sin_rot_angles = new double[len - 1];
        //flag = new int [2];
        //flag[0] = 0;
        //flag[1] = 0;

        for(j = 0; j < len-1;j++){
            cos_rot_angles[j] = A.cos_rot_angles[j];
            sin_rot_angles[j] = A.sin_rot_angles[j];
        
        }

        for(i = 0; i < _amountThreads; i++){
            _args[i].cos_rot_angles_str = cos_rot_angles;
            _args[i].sin_rot_angles_str = sin_rot_angles;
            _args[i].number_task = i;
            _args[i].number_thread = i;
            //_args[i].flag = flag;
            _args[i]._t_full = &_t_full;

        }

        
}

/*_____________________________________*/
Matrix::Matrix(int n, int __amountThreads){
        int i;

        _amountThreads = __amountThreads;
        len = n;
        
        _matrix = new double[len*len];
        _threads = new pthread_t[_amountThreads];  // массив потоков
        _args = new optionalStruct[_amountThreads]; // массив структур


        cos_rot_angles = new double[len - 1];
        sin_rot_angles = new double[len - 1];
        //flag = new int [1];
        //flag[0] = 0;
        //flag[1] = 0;

        for(i = 0; i < _amountThreads; i++){
            _args[i].cos_rot_angles_str = cos_rot_angles;
            _args[i].sin_rot_angles_str = sin_rot_angles;
            _args[i].number_task = i;
            _args[i].number_thread = i;
            //_args[i].flag = flag;
            _args[i]._t_full = &_t_full;
        }

}

/*_____________________________________*/
Matrix::Matrix(double *mass, int n, int __amountThreads){
        
        int i, j;
        _amountThreads = __amountThreads;
        len = n;

        _matrix = new double[len*len];
        for(i = 0; i < len; i++){
            for( j = 0; j < len; j++){
                _matrix[i*n+j] = mass[i*n+j];
            }   
        }
    
        _threads = new pthread_t[_amountThreads];  // массив потоков
        _args = new optionalStruct[_amountThreads]; // массив структур

        
        cos_rot_angles = new double[len - 1];
        sin_rot_angles = new double[len - 1];
        //flag = new int [1];
        //flag[0] = 0;
        //flag[1] = 0;

        for(i = 0; i < _amountThreads; i++){
            _args[i].cos_rot_angles_str = cos_rot_angles;
            _args[i].sin_rot_angles_str = sin_rot_angles;
            _args[i].number_task = 0;
            _args[i].number_thread = i;
            //_args[i].flag = flag;
            _args[i]._t_full = &_t_full;
        }
        

}

/*_____________________________________*/
Matrix::~Matrix(){

    
   // delete [] flag;
    delete [] cos_rot_angles;
    delete [] sin_rot_angles;
    delete [] _matrix;
    delete [] _threads;
    delete [] _args;

}

/*_____________________________________*/
int Matrix::InMatrix_formulas(int k){
    int i, j;
    switch (k)
    {
    case 1:
        for(i = 0; i < len; i++){
            for(j = 0; j < len; j++){
                _matrix[i*(len) +j] = len - max(i+1, j+1) + 1;
            }
        }
        return 0;
        break;
    case 2:
        for(i = 0; i < len; i++){
            for(j = 0; j < len; j++){
                _matrix[i*(len) +j] = max(i+1, j+1);
            }
        }
        return 0;
        break;
    case 3:
        for(i = 0; i < len; i++){
            for(j = 0; j < len; j++){
                _matrix[i*(len) +j] = abs(i-j);
            }
        }
        return 0;
        break;
    case 4:
        double key;
        for(i = 0; i < len; i++){
            for(j = 0; j < len; j++){
                key = 1.0/(i+j+1);
                _matrix[i*(len) +j] = key;
            }
        }
        return 0;
        break;
    
    default:
        return -1;
        break;
    }
}

/*_____________________________________*/
int Matrix::InMatrix_file(string nameOfFile){
    int i, j, kol = 0;
    int k;
    double x;
    FILE *f;
    if(!(f = fopen(nameOfFile.c_str(), "r"))){return -1;}
    
    while(fscanf(f, "%lf", &x) != -1){kol++;}
    if(kol != len*len){return -1;}
    
    fclose(f);

    if(!(f = fopen(nameOfFile.c_str(), "r"))){return -1;}

    for(i = 0; i < len; i++){
        for(j = 0; j < len; j++){
            k = fscanf(f, "%lf", &x);
            if(k == -1){
                fclose(f);
               return -1;
            }
            _matrix[i*(len) +j] = x;
        }
    }
    

    fclose(f);
    return 0;
}

/*_____________________________________*/
void Matrix::Print(int k){

    int i,j;
    for(i = 0; i < k; i++){
        for(j = 0; j < k; j++){
            printf(" %10.3e ", _matrix[i*(len) + j]);
        }
        cout << endl;
    }
}

static double total_time = 0.;
static pthread_mutex_t total_mutex = PTHREAD_MUTEX_INITIALIZER;

/*_____________________________________*/
void *process_function(void *pa){
    
    optionalStruct *temp = (optionalStruct*) pa;
    long int t;


    t = get_full_time();
    synchronize(temp->_amountThreads);
    parallel_Rotation_Method(temp->len, temp->copy_matrix, temp->inv_matrix, 
                            temp->number_thread, temp->_amountThreads,
                            temp->sin_rot_angles_str, temp->cos_rot_angles_str);
    synchronize(temp->_amountThreads);
    t = get_full_time() - t;
    
    *temp->_t_full = t;

    pthread_mutex_lock(&total_mutex);
    total_time += t;
    pthread_mutex_unlock(&total_mutex);


    
    return 0;
}

/*_____________________________________*/
int Matrix:: paralel_inverse_matrix(double *inv_matrix, double *copy_matrix){
    
    int i, j;
    double sum = 0.;
    double sub_const;
    
    /*_______________________________________________________*/
    /*The matrix is normalized*/
    
    for(i = 0; i < len; i++){
        for(j = 0; j < len; j++){
            sum = sum + (_matrix[i*(len) + j])*(_matrix[i*(len) + j]);
        }
    }
    sum = sqrt(sum);

    for(i = 0; i < len; i++){
        for(j = 0; j < len; j++){
            copy_matrix[i*len + j] = _matrix[i*len + j]/sum;
        }
    }
    
   

    /*_______________________________________________________*/
    /*A single matrix is created*/
    for(i = 0; i < len; i++){
        for(j = 0; j < len; j++){
            if(i == j){
                inv_matrix[i*len + j] = 1.;
            }else{
                inv_matrix[i*len + j] = 0.;
            }
        }
    }

    for(i = 0; i < _amountThreads; i ++){
        _args[i].copy_matrix = copy_matrix;
        _args[i].inv_matrix = inv_matrix;
        _args[i].len = len;
        _args[i]._amountThreads = _amountThreads;

    }

    

    /*_______________________________________________________*/
    
    /*The beginning of the method of Rotations reduction to a diagonal form*/
    printf("Start -> ");

    for(i = 0; i < _amountThreads; i++){
            if(pthread_create(_threads + i, 0, process_function, _args + i)){
                printf("Matrix error: cant create pthread_create\n");
                return -1;
            }
    }


    for(i = 0; i < _amountThreads; ++i){
        if(pthread_join (_threads[i], 0)){
            printf("Matrix error: cant wait pthread_create\n");
            return -2;
        }
    }

    printf(" end\n");
    





    /*_______________________________________________________*/
    /*Matrix normalization*/
    for(i = 0; i < len; i++){
        sub_const = 1./copy_matrix[i*len+i];
        for(j = 0; j < len; j++){
            inv_matrix[i*len+j] = (inv_matrix[i*len+j]*sub_const)/sum;
        }
    }



    return 0;

}

/*_____________________________________*/
double Restnormen(const Matrix &A, const Matrix &B){
    int i, j;
    int _len;
    double maxx = 0.;
    double sum = 0.;
    double res;

    if(A.len != B.len){
        return -1.;
    }

    _len = A.len;

    for(i = 0; i < _len; i++){
        for(j = 0; j < _len; j++){
            sum = A._matrix[i*_len + j]*B._matrix[j*_len+i] + sum;
        }

        res = sum - 1.0;
        if(res > maxx){
            maxx = res;
        }
        sum = 0.;

    }

    return maxx;
    
}

long int get_time (){
    struct rusage buf;
    getrusage (RUSAGE_SELF, &buf);
    return buf.ru_utime.tv_sec * 100 + buf.ru_utime.tv_usec / 10000;
}

long int get_full_time (){
    struct timeval buf;
    gettimeofday (&buf, 0);
    return buf.tv_sec * 100 + buf.tv_usec / 10000;
}

long int Matrix::get_elapsed(){
    return _t_full;
}



