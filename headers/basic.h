#ifndef CML_cml_basic_H  
#define CML_cml_basic_H


    #define NUM_LEN 8
    #define FLOAT_LEN 5

//    typedef enum {
//        INT,
//        FLOAT,
//        DOUBLE,
//        CHAR
//    } DataType;

	#ifndef CML_MATRIX_STRUCT_BODY
    	#define CML_MATRIX_STRUCT_BODY
        typedef struct __CML_MATRIX__ {
            int row_n;
            int col_n;
            double *m;
        } cml_Matrix_t;
    #endif




    void cml_basicFreeMatrix(cml_Matrix_t **Mp);

    int cml_basicGetIndex(cml_Matrix_t *M, int row, int col);

    void cml_basicShowMatrix(cml_Matrix_t *M);

    void cml_basicFillMatrixFromArray(double *array, cml_Matrix_t *M);

    cml_Matrix_t* cml_basicArray2Matrix(double *array, int row, int col);

    cml_Matrix_t* cml_basicCopy(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicRange(int from, int to, int skip);

    cml_Matrix_t* cml_basicEmpty(int row, int col);

    cml_Matrix_t* cml_basicZeros(int row, int col);

    cml_Matrix_t* cml_basicOnes(int row, int col);

    cml_Matrix_t* cml_basicZerosLike(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicOnesLike(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicIdentity(int dim);

    cml_Matrix_t* cml_basicIdentityLike(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicTranspose(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicAdd(cml_Matrix_t *A, cml_Matrix_t *B);

    cml_Matrix_t* cml_basicMinus(cml_Matrix_t *A, cml_Matrix_t *B);

    cml_Matrix_t* cml_basicNumProd(double a, cml_Matrix_t *A);

    cml_Matrix_t* cml_basicDot(cml_Matrix_t *A, cml_Matrix_t *B);

    cml_Matrix_t* __cml_basicSubMatrix__(cml_Matrix_t* M, int row, int col);

    double cml_basicDeterminant(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicInverse(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicDiag(double *elem, int dim);

    cml_Matrix_t* cml_basic_3d_cross(cml_Matrix_t *A, cml_Matrix_t *B);

    cml_Matrix_t* cml_basicSlice(cml_Matrix_t *M, cml_Matrix_t *RowSlice, cml_Matrix_t *ColSlice);

    cml_Matrix_t* cml_basicRSlice(cml_Matrix_t *M, cml_Matrix_t *RowSlice);

    cml_Matrix_t* cml_basicCSlice(cml_Matrix_t *M, cml_Matrix_t *ColSlice);

    cml_Matrix_t* cml_basicSum(cml_Matrix_t *M, int axis);

    cml_Matrix_t* cml_basicMean(cml_Matrix_t *M, int axis);

    cml_Matrix_t* cml_basicMatPow(cml_Matrix_t *M, int p);

    cml_Matrix_t* cml_basicElemPow(cml_Matrix_t *M, double p);

    cml_Matrix_t* cml_basicVecSort(cml_Matrix_t *Vec);

    cml_Matrix_t* cml_basicMatSort(cml_Matrix_t *M, int axis, int index);

    cml_Matrix_t* cml_basicReverse(cml_Matrix_t *Vec);

    double cml_basicVecMax(cml_Matrix_t *Vec);

    double cml_basicVecMin(cml_Matrix_t *Vec);

    cml_Matrix_t* cml_basicMatMax(cml_Matrix_t *M, int axis);

    cml_Matrix_t* cml_basicMatMin(cml_Matrix_t *M, int axis);

    int cml_basicApply(cml_Matrix_t *From, cml_Matrix_t *To, cml_Matrix_t *RowPos, cml_Matrix_t *ColPos);

    cml_Matrix_t* cml_basicShape(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicReshape(cml_Matrix_t *M, cml_Matrix_t *Shape);

    cml_Matrix_t* cml_basicFlatten(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicConcatenate(cml_Matrix_t *A, cml_Matrix_t *B, int axis);


#endif
