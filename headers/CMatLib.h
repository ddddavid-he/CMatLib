
#ifndef CML_C_MAT_LIB_H
#define CML_C_MAT_LIB_H 

    typedef struct __MATRIX_OBJ__ Mat;
    typedef struct __CML_MATRIX__ Matrix;

    #ifndef CML_BASIC_H 
    #ifndef CML_MATRIX_STRUCT_BODY
    #define CML_MATRIX_STRUCT_BODY
        typedef struct __CML_MATRIX__ {
            int row_n;
            int col_n;
            double *m;
        } Matrix;
    #endif
    #endif

    #define ERROR 0
    #define SUCCEED 1
    #define WARN -1

    #define TRUE 1
    #define FALSE 0

    typedef char bool;



    int cml_freeMatrix(Matrix **Mp);

    int cml_getIndex(Matrix *M, int row, int col);

    int cml_showMatrix(Matrix *M);

    int cml_fillMatrixFromArray(double *array, int row, int col, Matrix *M);

    Matrix* cml_array2Matrix(double *array, int row, int col);

    Matrix* cml_copy(Matrix *M);

    Matrix* cml_range(int from, int to, int step);

    Matrix* cml_zeros(int row, int col);

    Matrix* cml_ones(int row, int col);

    Matrix* cml_zerosLike(Matrix *M);

    Matrix* cml_onesLike(Matrix *M);

    Matrix* cml_identity(int dim);

    Matrix* cml_transpose(Matrix *M);

    Matrix* cml_add(Matrix *A, Matrix *B);

    Matrix* cml_minus(Matrix *A, Matrix *B);

    Matrix* cml_numProd(double a, Matrix *M);

    Matrix* cml_dot(Matrix *A, Matrix *B);

    double cml_det(Matrix *M);

    Matrix* cml_inv(Matrix *M);

    Matrix* cml_diag(double *elem, int dim);

    Matrix* cml_cross_3d(Matrix *A, Matrix *B);

    Matrix* cml_slice(Matrix *M, Matrix *Row_Idx, Matrix *Col_Idx);

    Matrix* cml_rowSlice(Matrix *M, Matrix *Row_Idx);

    Matrix* cml_colSlice(Matrix *M, Matrix *Col_Idx);

    Matrix* cml_sum(Matrix *M, int axis);

    Matrix* cml_mean(Matrix *, int axis);

    Matrix* cml_power(Matrix *M, int p);

    Matrix* cml_power_elementWise(Matrix *M, double p);

    Matrix* cml_sortVec_m(Matrix *Vec);

    Matrix* cml_sort_m(Matrix *M, int axis, int index);

    Matrix* cml_reverse_m(Matrix *Vec);

    double cml_vecMax_m(Matrix *Vec);

    double cml_vecMin_m(Matrix *Vec);

    Matrix* cml_max_m(Matrix *M, int axis);

    Matrix* cml_min_m(Matrix *M, int axis);

    int cml_apply_m(Matrix *From, Matrix *To, Matrix *RowPosition, Matrix *ColPosition);

    Matrix* cml_shape_m(Matrix *M);

    Matrix* cml_reshape_m(Matrix *M, Matrix *Shape);

    Matrix* cml_reshape2_m(Matrix *M, int row, int col);

    Matrix* cml_flatten_m(Matrix *M);

    // TODO:
    Matrix* cml_concatenate(Matrix *A, Matrix *B, int axis);

#endif
