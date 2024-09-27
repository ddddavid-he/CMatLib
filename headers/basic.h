#ifndef CML_BASIC_H  
#define CML_BASIC_H


    #define NUM_LEN 8
    #define FLOAT_LEN 5

//    typedef enum {
//        INT,
//        FLOAT,
//        DOUBLE,
//        CHAR
//    } DataType;

    typedef struct __MATRIX__ {
        int row_n;
        int col_n;
        double *m;
        DataType dtype;
    } Matrix;




    void basicFreeMatrix(Matrix **Mp);

    int giveIndex(Matrix *M, int row, int col);

    void basicShowMatrix(Matrix *M);

    void basicFillMatrixFromArray(double *array, int row, int col, Matrix *M);

    Matrix* basicArray2Matrix(double *array, int row, int col);

    Matrix* basicCopy(Matrix *M);

    Matrix* basicRange(int from, int to, int skip);

    Matrix* basicEmpty(int row, int col);

    Matrix* basicZeros(int row, int col);

    Matrix* basicOnes(int row, int col);

    Matrix* basicZerosLike(Matrix *M);

    Matrix* basicOnesLike(Matrix *M);

    Matrix* basicIdentity(int dim);

    Matrix* basicIdentityLike(Matrix *M);

    Matrix* basicTranspose(Matrix *M);

    Matrix* basicAdd(Matrix *A, Matrix *B);

    Matrix* basicMinus(Matrix *A, Matrix *B);

    Matrix* basicNumProd(double a, Matrix *A);

    Matrix* basicDot(Matrix *A, Matrix *B);

    Matrix* __subMatrix__(Matrix* M, int row, int col);

    double basicDeterminant(Matrix *M);

    Matrix* basicInverse(Matrix *M);

    Matrix* basicDiag(double *elem, int dim);

    Matrix* basic_3d_cross(Matrix *A, Matrix *B);

    Matrix* basicSlice(Matrix *M, Matrix *RowSlice, Matrix *ColSlice);

    Matrix* basicRSlice(Matrix *M, Matrix *RowSlice);

    Matrix* basicCSlice(Matrix *M, Matrix *ColSlice);

    Matrix* basicSum(Matrix *M, int axis);

    Matrix* basicMean(Matrix *M, int axis);

    Matrix* basicMatPow(Matrix *M, int p);

    Matrix* basicElemPow(Matrix *M, double p);

    Matrix* basicVecSort(Matrix *Vec);

    Matrix* basicMatSort(Matrix *M, int axis, int index);

    Matrix* basicReverse(Matrix *Vec);

    double basicVecMax(Matrix *Vec);

    double basicVecMin(Matrix *Vec);

    Matrix* basicMatMax(Matrix *M, int axis);

    Matrix* basicMatMin(Matrix *M, int axis);

    int basicApply(Matrix *From, Matrix *To, Matrix *RowPos, Matrix *ColPos);

    Matrix* basicShape(Matrix *M);

    Matrix* basicReshape(Matrix *M, Matrix *Shape);

    Matrix* basicFlatten(Matrix *M);

    Matrix* basicConcatenate(Matrix *A, Matrix *B, int axis);


#endif
