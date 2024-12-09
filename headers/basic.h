#ifndef CML_cml_basic_H  
#define CML_cml_basic_H

    // cml_showMatrix format definitions
    #define NUM_LEN 8
    #define FLOAT_LEN 5
    // end

	#ifndef CML_MATRIX_STRUCT_BODY
    	#define CML_MATRIX_STRUCT_BODY
        typedef struct __CML_MATRIX__ {
            int row_n;
            int col_n;
            double *m;
        } cml_Matrix_t;
    #endif





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
    void cml_basicFreeMatrix(cml_Matrix_t **Mp);


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
    int cml_basicGetIndex(cml_Matrix_t *M, int row, int col);


    /*
    Function to print a Matrix
    M: cml_Matrix_t*: Pointer of a Matrix
    return: void.
    Example:
    cml_Matrix_t *M;
    ...
    cml_basicShowMatrix(M);
    */
    void cml_basicShowMatrix(cml_Matrix_t *M);


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
    void cml_basicFillMatrixFromArray(double *array, cml_Matrix_t *M);


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
    cml_Matrix_t* cml_basicArray2Matrix(double *array, int row, int col);


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
    cml_Matrix_t* cml_basicCopy(cml_Matrix_t *M);


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
    cml_Matrix_t* cml_basicRange(int from, int to, int step);


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
    cml_Matrix_t* cml_basicEmpty(int row, int col);


    /*
    Create a Matrix of zeros in the shape of (row, col).
    row: int: row number of the Matrix.
    col: int: column number of the Matrix.
    return: cml_Matrix_t*: pointer of the Matrix.
    Example:
        cml_Matrix_t *M = cml_basicZeros(2, 1);
    // M == [ [0], [0] ]
    */
    cml_Matrix_t* cml_basicZeros(int row, int col);


    /*
    Create a Matrix of ones in the shape of (row, col);
    row: int: row number of the Matrix.
    col: int: column number of the Matrix.
    return: cml_Matrix_t*: pointer of the Matrix.
    Example:
        cml_Matrix_t *M = cml_basicOnes(2, 1);
    // M == [ [1], [1] ]
    */
    cml_Matrix_t* cml_basicOnes(int row, int col);


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
    cml_Matrix_t* cml_basicZerosLike(cml_Matrix_t *M);


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
    cml_Matrix_t* cml_basicOnesLike(cml_Matrix_t *M);


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
    cml_Matrix_t* cml_basicIdentity(int dim);


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
    cml_Matrix_t* cml_basicIdentityLike(cml_Matrix_t *M);


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
    cml_Matrix_t* cml_basicTranspose(cml_Matrix_t *M);


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
    cml_Matrix_t* cml_basicAdd(cml_Matrix_t *A, cml_Matrix_t *B);


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
    cml_Matrix_t* cml_basicMinus(cml_Matrix_t *A, cml_Matrix_t *B);


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
    cml_Matrix_t* cml_basicNumProd(double a, cml_Matrix_t *A);


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
    cml_Matrix_t* cml_basicDot(cml_Matrix_t *A, cml_Matrix_t *B);


    /*
    Function of getting the sub-matrix of a (row, col) of a Matrix.
    In linear algebra, A_{row,col}.
    M: cml_Matrix_t*: the origin Matrix.
    row, col: int: row index and col index of the sub-matrix (starting from 0)
    returning: cml_Matrix_t*: pointer of the returning sub-matrix
    Example:
        cml_Matrix_t *M, *sub_M;
        ...  // M == [ [1, 2, 3], [4, 5, 6], [7, 8, 9] ]
        sub_M = cml_basicSubMatrix(M, 1, 1);
    // sub_M == [ [1, 3], [7, 9] ]
    */
    cml_Matrix_t* cml_basicSubMatrix(cml_Matrix_t* M, int row, int col);


    /*
    Function of calculating determinant of a Matrix. det(M) in linear algebra.
    M: cml_Matrix_t*: pointer of the Matrix.
    return: value of determinant in precision of `double`.
    Example:
        cml_Matrix_t *M;
        ... // M == [ [1, 2], [3, 4] ]
        double ret = cml_basicDeterminant(M);
    // ret == 1*4 - 2*3 = -2
    */
    double cml_basicDeterminant(cml_Matrix_t *M);


    /*
    Calculating inverse matrix of a Matrix.
    In linear algebra, writes M^{-1}. M^{-1} * M == I;
    Example:
        cml_Matrix_t *A, *A_inv;
        ... // A == [ [1, 2], [3, 4] ]
        A_inv = cml_basicInverse(A);
    // A_inv == [ [-2, 1], [1.5, -0.5] ]
    */
    cml_Matrix_t* cml_basicInverse(cml_Matrix_t *M);


    /*
    Create a diagonal matrix from double pointer.
    Example:
        cml_Matrix_t *M;
        double *array = {1, 2, 3};
        cml_Matrix_t *A;
        A = cml_basicDiag(array, 3);
    // A == [ [1, 0, 0], [0, 2, 0], [0, 0, 3] ]
    */
    cml_Matrix_t* cml_basicDiag(double *elem, int dim);

    /*
    Perform cross product of two 3D vectors.
    Example:
        cml_Matrix_t *a;
        cml_Matrix_t *b;
        cml_Matrix_t *c;
        double *values_of_a = {1, 0, 0};
        double *values_of_b = {0, 1, 0};
        a = cml_basicArray2Matrix(values_of_a, 1, 3);
        b = cml_basicArray2Matrix(values_of_b, 1, 3);
        c = cml_basicCross(a, b);
    // c == [ [0, 0, 1] ]
    */
    cml_Matrix_t* cml_basicCross(cml_Matrix_t *A, cml_Matrix_t *B);

    /*
    Create a slice of the given matrix by row indices and column indices.
    Example:
        cml_Matrix_t *A, *Row, *Col, *Slice;
        double *A_val = {0, 1, 2, 3};
        double *r_val = {0, 1};
        double *c_val = {1, 2};
        A = cml_basicDiag(A_val, 4);
        Row = cml_basicArray2Matrix(r_val, 1, 2);
        Col = cml_basicArray2Matrix(c_val, 1, 2);
        Slice = cml_basicSlice(A, Row, Col);
    // Slice == [ [0, 0], [1, 0] ]
    */
    cml_Matrix_t* cml_basicSlice(cml_Matrix_t *M, cml_Matrix_t *RowIdx, cml_Matrix_t *ColIdx);

    /*
    Create a row slice of the given matrix by row indices.
    Example:
        cml_Matrix_t *A, *Row, *Slice;
        double *A_val = {0, 1, 2, 3};
        double *r_val = {0, 1};
        A = cml_basicDiag(A_val, 4);
        Row = cml_basicArray2Matrix(r_val, 1, 2);
        Slice = cml_basicRowSlice(A, Row);
    // Slice == [ [0, 0, 0, 0], [0, 1, 0, 0] ]
    */
    cml_Matrix_t* cml_basicRowSlice(cml_Matrix_t *M, cml_Matrix_t *RowSlice);

    /*
    Create a column slice of the given matrix by column indices.
    Example:
        cml_Matrix_t *A, *Col, *Slice;
        double *A_val = {0, 1, 2, 3};
        double *c_val = {3, 2};
        A = cml_basicDiag(A_val, 4);
        Col = cml_basicArray2Matrix(c_val, 1, 2);
        Slice = cml_basicColSlice(A, Col);
    // Slice == [ [0, 0], [0, 0], [0, 2], [3, 0] ]
    */
    cml_Matrix_t* cml_basicColSlice(cml_Matrix_t *M, cml_Matrix_t *ColSlice);

    /*
    Sum of the elements along `axis`.
        axis==-1: sum of all elements;
        axis==0 : sum of each row;
        axis==1 : sum of each column;
    Example:
        array *values = {1,2,3,4};
        cml_Matrix_t *M = cml_basicArray2Matrix(values, 2, 2);
        int s_a, s_b, s_c;
        s_a = cml_basicSum(M, -1);
        s_b = cml_basicSum(M, 0);
        s_c = cml_basicSum(M, 1);
    // s_a == 10
    // s_b == [ [3], [7] ]
    // s_c == [ [4, 6] ];
    */
    cml_Matrix_t* cml_basicSum(cml_Matrix_t *M, int axis);

    /*
    Mean value of elements along `axis`.
    axis == -1: mean value of all elements;
    axis == 0 : mean value of each row;
    axis == 1 : mean value of each column;
    Example:
        double *values = {1,2,3,4};
        cml_Matrix_t *M = cml_basicArray2Matrix(values, 2, 2);
        int m_a, m_b, m_c;
        m_a = cml_basicSum(M, -1);
        m_b = cml_basicSum(M, 0);
        m_c = cml_basicSum(M, 1);
    // m_a == 2.5
    // m_b == [ [1.5], [3.5] ]
    // m_c == [ [2, 3] ];
    */
    cml_Matrix_t* cml_basicMean(cml_Matrix_t *M, int axis);

    /*
    Matrix of M^p (dot product of p Ms).
    Example:
        double *array = {1, 2};
        cml_Matrix_t *M = cml_basicDiag(array, 2);
        cml_Matrix_t *res;
        res = cml_basicMatPow(M, 3);
    // res == [ [1, 0], [0, 8] ];
    */
    cml_Matrix_t* cml_basicMatPow(cml_Matrix_t *M, int p);

    /*
    Power operation for each element.
    Example:
        double *values = {0, 1, 2, 3};
        cml_Matrix_t *M = cml_basicArray2Matrix(values, 2, 2);
        cml_Matrix_t *res;
        res = cml_basicElemPow(M, 2);
    // res == [ [0, 1], [4, 9] ]
    */
    cml_Matrix_t* cml_basicElemPow(cml_Matrix_t *M, double p);

    /*
    Sort values in a vector. (Ascending, returning a copy)
    Example:
        cml_Matrix_t *M, *sorted;
        double *array = {2,4,3,1};
        M = cml_basicArray2Matrix(array, 4, 1);
        sorted = cml_basicVecSort(M);
    // M == [ [2, 4, 3, 1] ]
       sorted == [ [1, 2, 3, 4] ]
    */
    cml_Matrix_t* cml_basicVecSort(cml_Matrix_t *Vec);

    /*
    Sort rows or columns in the order of values from a vector specified by `axis` and `index`.
    axis == 0: sort by row number (`index`)
    axis == 1: sort by col number (`index`)
    Example:
        M == [[4,5,6], [2,1,0], [1,3,3]]
        axis=0, index=1: sorted = [[6,5,4], [0,1,2], [3,3,1]]
        axis=1, index=0: sorted = [[1,3,3], [2,1,0], [4,5,6]]
    */
    cml_Matrix_t* cml_basicMatSort(cml_Matrix_t *M, int axis, int index);

    /*
    Returns a copy of a vector in its reversed order.
    Example:
        // M == [[ 0, 1, 2, 3 ]]
        cml_Matrix_t *M_r = cml_basicReverse(M);
    // M_r == [[ 3, 2, 1, 0 ]]
    */
    cml_Matrix_t* cml_basicReverse(cml_Matrix_t *Vec);

    /*
    Returns the maximum value in a vector or matrix.
    Example:
        // M == [ [1, 2, 3], [4, 5, 6] ]
        double max = cml_basicVecMax(M);
    // max == 1
    */
    double cml_basicVecMax(cml_Matrix_t *Vec);

    /*
    Returns the minimum value in a vector or matrix.

    Example:
        // M == [ [1, 2, 3], [4, 5, 6] ]
        double max = cml_basicVecMax(M);
    // max == 1
    */
    double cml_basicVecMin(cml_Matrix_t *Vec);

    /*
    Returns maximum values of rows or columns of a matrix.
    M: cml_Matrix_t*: the input matrix.
    axis: int: 0 or 1, calculate maximum value along axis `axis`.
    axis == 0: returns a vector of (row, 1) of maximum values of each row.
    axis == 1: returns a vector of (1, col) of maximum values of each column.
    Example:
        // M == [ [1, 2, 3], [4, 5, 6] ]
        cml_Matrix_t* ret0 = cml_basicMatMax(M, 0);
        cml_Matrix_t* ret1 = cml_basicMatMax(M, 1);
    // ret0 == [ [3], [6] ]
       ret1 == [ [4, 5, 6] ]
    */
    cml_Matrix_t* cml_basicMatMax(cml_Matrix_t *M, int axis);

    /*
    Returns minimum values of row or columns of a matrix.
    M: cml_Matrix_t*: the input matrix.
    axis: int: 0 or 1, calculate minimum value along axis `axis`.
    axis == 0: returns a vector of (row, 1) of minimum values of each row.
    axis == 1: returns a vector of (1, col) of minimum values of each column.
    Example:
        // M == [ [1, 2, 3], [4, 5, 6] ]
        cml_Matrix_t* ret0 = cml_basicMatMin(M, 0);
        cml_Matrix_t* ret1 = cml_basicMatMin(M, 1);
    // ret0 == [ [1], [4] ]
       ret1 == [ [1, 2, 3] ]
    */
    cml_Matrix_t* cml_basicMatMin(cml_Matrix_t *M, int axis);

    int cml_basicApply(cml_Matrix_t *From, cml_Matrix_t *To, cml_Matrix_t *RowPos, cml_Matrix_t *ColPos);

    cml_Matrix_t* cml_basicShape(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicReshape(cml_Matrix_t *M, cml_Matrix_t *Shape);

    cml_Matrix_t* cml_basicFlatten(cml_Matrix_t *M);

    cml_Matrix_t* cml_basicConcatenate(cml_Matrix_t *A, cml_Matrix_t *B, int axis);


#endif
