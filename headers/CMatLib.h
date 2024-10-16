
#ifndef CML_C_MAT_LIB_H
#define CML_C_MAT_LIB_H 

    typedef struct __MATRIX_OBJ__ Mat;
//    typedef struct __CML_MATRIX__ Matrix;

    #ifndef CML_BASIC_H 
        #ifndef CML_MATRIX_STRUCT_BODY
        #define CML_MATRIX_STRUCT_BODY
            typedef struct __CML_MATRIX__ {
                int row_n;
                int col_n;
                double *m;
            } cml_Matrix_t;
        #endif
    #endif

    #define ERROR 0
    #define SUCCEED 1
    #define WARN -1

    #define TRUE 1
    #define FALSE 0

    typedef char bool;



    int cml_freeMatrix(cml_Matrix_t **Mp);

    int cml_getIndex(cml_Matrix_t *M, int row, int col);

    int cml_showMatrix(cml_Matrix_t *M);

    int cml_fillMatrixFromArray(double *array, int row, int col, cml_Matrix_t *M);

    cml_Matrix_t* cml_array2Matrix(double *array, int row, int col);

    cml_Matrix_t* cml_copy(cml_Matrix_t *M);

    cml_Matrix_t* cml_range(int from, int to, int step);

    cml_Matrix_t* cml_zeros(int row, int col);

    cml_Matrix_t* cml_ones(int row, int col);

    cml_Matrix_t* cml_zerosLike(cml_Matrix_t *M);

    cml_Matrix_t* cml_onesLike(cml_Matrix_t *M);

    cml_Matrix_t* cml_identity(int dim);

    cml_Matrix_t* cml_transpose(cml_Matrix_t *M);

    cml_Matrix_t* cml_add(cml_Matrix_t *A, cml_Matrix_t *B);

    cml_Matrix_t* cml_minus(cml_Matrix_t *A, cml_Matrix_t *B);

    cml_Matrix_t* cml_numProd(double a, cml_Matrix_t *M);

    cml_Matrix_t* cml_dot(cml_Matrix_t *A, cml_Matrix_t *B);

    double cml_det(cml_Matrix_t *M);

    cml_Matrix_t* cml_inv(cml_Matrix_t *M);

    cml_Matrix_t* cml_diag(double *elem, int dim);

    cml_Matrix_t* cml_cross_3d(cml_Matrix_t *A, cml_Matrix_t *B);

    cml_Matrix_t* cml_slice(cml_Matrix_t *M, cml_Matrix_t *Row_Idx, cml_Matrix_t *Col_Idx);

    cml_Matrix_t* cml_rowSlice(cml_Matrix_t *M, cml_Matrix_t *Row_Idx);

    cml_Matrix_t* cml_colSlice(cml_Matrix_t *M, cml_Matrix_t *Col_Idx);

    cml_Matrix_t* cml_sum(cml_Matrix_t *M, int axis);

    cml_Matrix_t* cml_mean(cml_Matrix_t *, int axis);

    cml_Matrix_t* cml_power(cml_Matrix_t *M, int p);

    cml_Matrix_t* cml_power_elementWise(cml_Matrix_t *M, double p);

    cml_Matrix_t* cml_sortVec(cml_Matrix_t *Vec);

    cml_Matrix_t* cml_sort(cml_Matrix_t *M, int axis, int index);

    cml_Matrix_t* cml_reverse(cml_Matrix_t *Vec);

    double cml_vecMax(cml_Matrix_t *Vec);

    double cml_vecMin(cml_Matrix_t *Vec);

    cml_Matrix_t* cml_max(cml_Matrix_t *M, int axis);

    cml_Matrix_t* cml_min(cml_Matrix_t *M, int axis);

    int cml_apply(cml_Matrix_t *From, cml_Matrix_t *To, cml_Matrix_t *RowPosition, cml_Matrix_t *ColPosition);

    cml_Matrix_t* cml_shape(cml_Matrix_t *M);

    cml_Matrix_t* cml_reshape(cml_Matrix_t *M, cml_Matrix_t *Shape);

    cml_Matrix_t* cml_reshape2(cml_Matrix_t *M, int row, int col);

    cml_Matrix_t* cml_flatten(cml_Matrix_t *M);

    cml_Matrix_t* cml_concatenate(cml_Matrix_t *A, cml_Matrix_t *B, int axis);

#endif
