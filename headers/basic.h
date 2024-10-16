#ifndef CML_cml_basic_H  
#define CML_cml_basic_H


    #define NUM_LEN 8
    #define FLOAT_LEN 5

	#ifndef CML_MATRIX_STRUCT_BODY
    	#define CML_MATRIX_STRUCT_BODY
        typedef struct __CML_MATRIX__ {
            int row_n;
            int col_n;
            double *m;
        } cml_Matrix_t;
    #endif




    void cml_basicFreeMatrix(cml_Matrix_t **Mp);
    /*
    Free up the space allocated to a pointer of Matrix, Matrix->m included.
    Mp: cml_Matrix_t**: Address of a cml_Matrix_t Pointer
    return: void
    Example:
     cml_Matrix_t *M;
     ...
     cml_basicFreeMatrix(&M);
    // M->m == NULL, M == NULL
    Notice: Mp should be the `address` of a `cml_Matrix_t Pointer`.
    */

    int cml_basicGetIndex(cml_Matrix_t *M, int row, int col);
    /*
    Convert 2D-index of (row, col) into 1D-index of M->m.
    M: cml_Matrix_t*: Which Matrix you want to get index of.
    i: int: row index.
    j: int: column index.
    return: int: 1D-index of M->m.
    Example:
        cml_Matrix_t M;
        M.row_n = 2;
        M.col_n = 3;
        int idx = 0, value=0;
        for(int i=0; i<M.row_n; i++) {
            for(int j=0; j<M.row_n; j++) {
                idx = cml_basicGetIndex(&M, i, j);
                M.m[idx] = value++
            }
        }
    // M == [ [0, 1, 2],
    //        [3, 4, 5] ]
    */

    void cml_basicShowMatrix(cml_Matrix_t *M);
    /*
    Function to print a Matrix
    M: cml_Matrix_t*: Pointer of a Matrix
    return: void.
    Example:
    cml_Matrix_t *M;
    ...
    cml_basicShowMatrix(M);
    */

    void cml_basicFillMatrixFromArray(double *array, cml_Matrix_t *M);
    /*
    Fill values from an array into a Matrix.
    Shape of the Matrix has to be defined and Matrix->m has to be allocated
    before passing Matrix to this function.
    array: double*: array of values to be filled.
    M: cml_Matrix_t*: Matrix to be filled.
    return: void.
    Example:
      cml_Matrix_t M;
      M.row_n = 2;
      M.col_n = 2;
      M.m = (double*) calloc(M.row_n*M.col_n, sizeof(double));
      double array[4] = {1, 2, 3, 4};
      cml_basicFillMatrixFromArray(array, &M);
    // M == [ [1, 2],
    //        [3, 4] ]
    */

    cml_Matrix_t* cml_basicArray2Matrix(double *array, int row, int col);
    /*
    Create a Matrix from array with specified (row, col).
    array: double*: array of values.
    row: int: row number of the Matrix.
    col: int: col number of the Matrix.
    return: cml_Matrix_t*: pointer of the Matrix.
    Example:
        cml_Matrix_t *M;
        double array[4] = {1, 2, 3, 4};
        M = cml_basicArray2Matrix(array, 2, 2);
    // M == [ [1, 2],
    //        [3, 4] ]
    */

    cml_Matrix_t* cml_basicCopy(cml_Matrix_t *M);
    /*
    Copy a Matrix to another.
    M: cml_Matrix_t*: source Matrix.
    return: cml_Matrix_t*: target Matrix.
    Example:
        cml_Matrix_t *old, *new;
        old = ...;
        ...
        new = cml_basicCopy(old);
    // new == old
    Notice: returning Matrix is re-allocated Matrix instead of a reference.
    */

    cml_Matrix_t* cml_basicRange(int from, int to, int step);
    /*
    Create a Matrix of arithmetic progression.
    from: int: start number of the series.
    to: int: where the series ends.
    step: int: step size.
    Example:
        cml_Matrix_t *series = cml_basicRange(1, 11, 2);
    // series == [[1, 3, 5, 7, 9]]
    Notice: `to` is not included while `from` is.
    */

    cml_Matrix_t* cml_basicEmpty(int row, int col);
    /*
    Create an empty matrix (allocated but not initialized).
    row: int: row number of the Matrix.
    col: int: col number of the Matrix.
    return: cml_Matrix_t*: pointer of the Matrix.
    Example:
        cml_Matrix_t *M = cml_basicEmpty(2, 2);
        M->m[0] = 1;
        M->m[0] = 2;
        M->m[0] = 3;
        M->m[0] = 4;
    // M == [ [1, 2], [3, 4] ]
    */

    cml_Matrix_t* cml_basicZeros(int row, int col);
    /*
    Create a Matrix of zeros in the shape of (row, col).
    row: int: row number of the Matrix.
    col: int: column number of the Matrix.
    return: cml_Matrix_t*: pointer of the Matrix.
    Example:
        cml_Matrix_t *M = cml_basicZeros(2, 1);
    // M == [ [0], [0] ]
    */

    cml_Matrix_t* cml_basicOnes(int row, int col);
    /*
    Create a Matrix of ones in the shape of (row, col);
    row: int: row number of the Matrix.
    col: int: column number of the Matrix.
    return: cml_Matrix_t*: pointer of the Matrix.
    Example:
        cml_Matrix_t *M = cml_basicOnes(2, 1);
    // M == [ [1], [1] ]
    */


    cml_Matrix_t* cml_basicZerosLike(cml_Matrix_t *M);
    /*
    Create a Matrix of zeros in the shape of another given Matrix.
    M: cml_Matrix_t*: the Matrix
    return: cml_Matrix_t*: pointer of the returning Matrix.
    Example:
        cml_Matrix_t A;
        A.row_n = 1;
        A.col_n = 2;
        cml_Matrix_t *M = cml_basicZerosLike(&A);
    // M == [ [0, 0] ]
    */

    cml_Matrix_t* cml_basicOnesLike(cml_Matrix_t *M);
    /*
    Create a Matrix of ones in the shape of another given Matrix.
    M: cml_Matrix_t*: the Matrix
    return: cml_Matrix_t*: pointer of the returning Matrix.
    Example:
        cml_Matrix_t A;
        A.row_n = 1;
        A.col_n = 2;
        cml_Matrix_t *M = cml_basicOnesLike(&A);
    // M == [ [1, 1] ]
    */

    cml_Matrix_t* cml_basicIdentity(int dim);
    /*
    Create an Identity Matrix of given dimension.
    dim: int: dimension of the Identity Matrix.
    return: cml_Matrix_t*: pointer of the returning Matrix.
    Example:
        cml_Matrix_t *I = cml_basicIdentity(3);
    // I == [ [1, 0, 0],
              [0, 1, 0],
              [0, 0, 1] ]
    */

    cml_Matrix_t* cml_basicIdentityLike(cml_Matrix_t *M);
    /*
    Create an Identity Matrix of same dimension as a given Matrix.
    M: cml_Matrix_t*: the given Matrix.
    return: cml_Matrix_t*: pointer of the returning Matrix.
    Example:
    cml_Matrix_t *M = cml_basicOnes(3, 3);
    cml_Matrix_t *I = cml_basicIdentityLike(M);
    // I == [ [1, 0, 0],
              [0, 1, 0],
              [0, 0, 1] ]
    // Notice: if M is not square-shaped. Dimension of the returning Matrix
    // equals to M->row_n.
    */

    cml_Matrix_t* cml_basicTranspose(cml_Matrix_t *M);
    /*
    Transpose the given Matrix.
    (The returning value is a copy of the original Matrix)
    M: cml_Matrix_t*: the Matrix to be transposed.
    return: cml_Matrix_t*: pointer of returning Matrix.
    Example:
        double array[4] = {1, 2, 3, 4};
        cml_Matrix_t *M = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *Mt = cml_basicTranspose(M);
    // M == [ [1, 2],
              [3, 4] ];
       Mt == [ [1, 3],
               [2, 4] ];
    */

    cml_Matrix_t* cml_basicAdd(cml_Matrix_t *A, cml_Matrix_t *B);
    /*
    Matrix add function. Performs Matrix_A + Matrix_B.
    Matrices are added element-wised.
    A, B: cml_Matrix*: Matrix to be added.
    return: cml_Matrix_t*: pointer of returning Matrix.
    Example:
        double array[4] = {1, 2, 3, 4};
        cml_Matrix_t *A = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *B = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *M = cml_basicAdd(A, B);
    // M == [ [2, 4], [6, 8] ];
    */

    cml_Matrix_t* cml_basicMinus(cml_Matrix_t *A, cml_Matrix_t *B);
    /*
    Matrix minus function. Performs Matrix_A - Matrix_B element-wise.
    A, B: cml_Matrix*: Matrix to be calculated.
    return: cml_Matrix_t*: pointer of returning Matrix.
    Example:
        double array[4] = {1, 2, 3, 4};
        cml_Matrix_t *A = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *B = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *M = cml_basicMinus(A, B);
    // M == [ [0, 0], [0, 0] ];
    */

    cml_Matrix_t* cml_basicNumProd(double a, cml_Matrix_t *A);
    /*
    Function of multiplying a Matrix by a number. Performs a * Matrix_A.
    a: double: multiplying number.
    A: cml_Matrix*: Matrix to be multiplied.
    return: cml_Matrix_t*: pointer of returning Matrix.
    Example:
        double array[4] = {1, 2, 3, 4};
        cml_Matrix_t *A = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *M = cml_basicNumProd(2, A);
    // M == [ [2, 4], [6, 8] ];
    */

    cml_Matrix_t* cml_basicDot(cml_Matrix_t *A, cml_Matrix_t *B);
    /*
    Function of performing dot product of two matrices. Performs Matrix_A * Matrix_B.
    A, B: cml_Matrix*: Matrices to be multiplied.
    return: cml_Matrix_t*: pointer of returning Matrix.
    Example:
        double array[4] = {1, 2, 3, 4};
        cml_Matrix_t *A = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *B = cml_basicArray2Matrix(array, 2, 2);
        cml_Matrix_t *M = cml_basicDot(A, B);
    // M == [ [7, 10], [15, 22] ];
    */

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
