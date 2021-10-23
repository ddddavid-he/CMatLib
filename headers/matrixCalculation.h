
#ifndef MATRIX_CALCULATION_H 
#define MATRIX_CALCULATION_H


    typedef struct __MATRIX_OBJ__ Mat;
    typedef struct __MATRIX__ Matrix;

    #ifndef BASIC_CALCULATION_H 
    #ifndef MATRIX_STRUCT_BODY
    #define MATRIX_STRUCT_BODY
        typedef struct __MATRIX__ {
            int row_n;
            int col_n;
            double *m;
        } Matrix;
    #endif
    #endif



    int freeMatrix(Matrix **Mp);

    int getIndex(Matrix *M, int row, int col);

    int showMatrix(Matrix *M);

    int fillMatrixFromArray(double *array, int row, int col, Matrix *M);

    Matrix* array2Matrix(double *array, int row, int col);

    Matrix* copy_m(Matrix *M);

    Matrix* range(int from, int to, int skip);

    Matrix* zeros_m(int row, int col);

    Matrix* ones_m(int row, int col);

    Matrix* zerosLike_m(Matrix *M);

    Matrix* onesLike_m(Matrix *M);

    Matrix* identity_m(int dim);

    Matrix* identityLike_m(Matrix *M);

    Matrix* trans_m(Matrix *M);

    Matrix* add_m(Matrix *A, Matrix *B);

    Matrix* minus_m(Matrix *A, Matrix *B);

    Matrix* nProd_m(double a, Matrix *A);

    Matrix* dot_m(Matrix *A, Matrix *B);

    double det_m(Matrix *M);

    Matrix* inv_m(Matrix *M);

    Matrix* diag_m(double *elem, int dim);

    Matrix* cross_3d_m(Matrix *A, Matrix *B);

    Matrix* slice_m(Matrix *M, Matrix *Row_Slice, Matrix *Col_Slice);

    Matrix* rowSlice_m(Matrix *M, Matrix *Row_Slice); 

    Matrix* colSlice_m(Matrix *M, Matrix *Col_Slice); 

    Matrix* sum_m(Matrix *M, int axis);

    Matrix* mean_m(Matrix *, int axis);

    Matrix* matPow_m(Matrix *M, int p);

    Matrix* elemPow_m(Matrix *M, double p);

    Matrix* sortVec_m(Matrix *Vec);

    Matrix* sort_m(Matrix *M, int axis, int index);

    Matrix* reverse_m(Matrix *Vec);

    double vecMax_m(Matrix *Vec);

    double vecMin_m(Matrix *Vec);

    Matrix* max_m(Matrix *M, int axis);

    Matrix* min_m(Matrix *M, int axis);

    int apply_m(Matrix *From, Matrix *To, Matrix *RowPosition, Matrix *ColPosition);

    Matrix* shape_m(Matrix *M);

    Matrix* reshape_m(Matrix *M, Matrix *Shape);

    Matrix* reshape2_m(Matrix *M, int row, int col);

    Matrix* flatten_m(Matrix *M);

#endif
