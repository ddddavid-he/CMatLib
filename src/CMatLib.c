#ifndef CML_C_MAT_LIB
#define CML_C_MAT_LIB

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "basic.c"

#include "CMatLib.h"


// Error function and Warning function
void __cml_logger(
        const char *func, const char *file, int line, int type,
        const char *message, ...
) {
    va_list args;
    va_start(args, message);
    if(type==0){
        const char *type_name = "Error";
    } else {
        const char *type_name = "Warning";
    }
    fprintf(stderr, "%s in func %s: (File: %s, Line: %d) ", type_name, func, file, line);
    vfprintf(stderr, message, args);
    fprintf(stderr, "\n");
    va_end(args);
    if(type==0){
        exit(-1);
    }
}

#define cml_error(message, ...) __cml_logger(__func__, __FILE__, __LINE__, 0, message)
#define cml_warning(message, ...) __cml_logger(__func__, __FILE__, __LINE__, 1, message)





int cml_freeMatrix(Matrix **Mp) {
    // free up pointer of a Matrix
    // CAUTION: 
    // Mp have to be Pointer of Matrix Pointer.
    // Usage:
    // cml_freeMatrix(&pointer_of_matrix)
    if(Mp==NULL) {
        cml_error("pointer of Matrix pointer Mp is NULL.");
        return ERROR;
    }
    if( (*Mp)==NULL ) {
        cml_warning("matrix pointer is already NULL.");
        return WARN;
    }

    basicFreeMatrix(Mp);
    return SUCCEED;
}


Matrix* cml_empty(uint row, uint col) {
    Matrix* M = (Matrix*) malloc(sizeof(Matrix));
    M->m = calloc(row*col, sizeof(double));
    return M;
}


int cml_fillMatrixFromArray(double *array, int row, int col, Matrix *M) {
    // fill values from array to matrix M
    // notice: Matrix->m should be allocated before this function
    if( array==NULL ) {
        cml_error("array is NULL.");
        return ERROR;
    }
    if( row<=0 || col<=0 ) {
        cml_error("row and col should be positive.");
        return ERROR;
    }
    if( M==NULL ) {
        cml_error("Matrix is NULL.");
        return ERROR;
    }
    if( M->m == NULL ) {
        cml_error("Matrix->m is not allocated yet.");
    }

    basicFillMatrixFromArray(array, row, col, M);
    return SUCCEED;
}


Matrix* cml_arrayToMatrix(double *array, int row, int col) {
    // create a Matrix from values of array
    Matrix *M = (Matrix*) malloc(sizeof(Matrix));
    M-> m = (double*) calloc(row*col, sizeof(double));
    cml_fillMatrixFromArray(array, row, col, M);
    return M;
}


Matrix* __cml_reformIndex__(int RowCol, Matrix *Indices) {
    // enable index to be in (-l, l-1)
    Matrix *idx = NULL;
    if(Indices->row_n!=1 && Indices->col_n==1){
        idx = basicTranspose(Indices);
    }else{
        idx = basicCopy(Indices);
    }

    for(int i=0;i<idx->row_n;i++){
        for(int j=0;j<Indices->col_n;j++){
            int index = giveIndex(in, i, j);
            if(idx->m[index] < 0){
                if(idx->m[index] < -RowCol){
                    cml_error("index %d out of range %d.", (int)idx->m[index], RowCol);
                    return NULL;
                }
                idx->m[index] += RowCol;
            }
        }
    }
    return idx;
}


int cml_getIndex(Matrix* M, int row, int col) {
    // transform row index and col index into 1D index
    // usage:
    // M->m[cml_getIndex(M, 2, 2)]
    if(M==NULL){
        cml_error("Matrix is NULL.");
    }
    if(M->m==NULL){
        cml_warning("Matrix->m is NULL.");
    }
    if(i >= M->row_n && i < -M->row_n){
        cml_error("getIndex: row index i out of range. Row=%d\n", row);
        return NULL
    }else if(j >= M->col_n && j < -M->col_n){
        cml_error("getIndex: column index j out of range. Col=%d\n", col);
        return NULL
    }

    if(row<0){
        row += M->row_n;
    }
    if(col<0){
        col += M->col_n;
    }
    return giveIndex(M, row, col);
}


int cml_showMatrix(Matrix *M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return ERROR;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL");
        return ERROR;
    }
    basicShowMatrix(M);
    return SUCCEED;
}


Matrix* cml_copy(Matrix *M) {
    if( M==NULL || M->m==NULL ){
        cml_error("Matrix is NULL.");
        return NULL;
    }

    return basicCopy(M);
}


Matrix* cml_range(int from, int to, int step) {
    if(from==to){
        if(step==0){
            return basicNumProd(from, basicOnes(1, 1));
        }else{
            cml_error("Empty range.");
            return NULL;
        }
    }else if(from<to){
        if(step<=0){
            cml_error("Emtpy range.");
            return NULL;
        }
        return basicRange(from, to, step);
    }else{ 
        // from > to
        if(step>=0){
            cml_error("Empty range.");
            return NULL;
        }
        return basicRange(from, to, step);
    }
}


Matrix* cml_zeros(int row, int col) {
    if(row*col<=0){
        cml_error("Cannot create a Matrix with non-positive row or column number.");
        return NULL;
    }
    return basicZeros(row, col);
}


Matrix* cml_ones(int row, int col) {
    if(row*col<=0){
        cml_error("Cannot create a Matrix with non-positive row or column number.");
        return NULL;
    }
    return basicOnes(row, col);
}


Matrix* cml_zerosLike(Matrix* M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    return cml_zeros(M->row_n, M->col_n);
}


Matrix* cml_onesLike(Matrix* M) {
    if(M==NULL) {
        cml_error("Matrix is NULL.");
        return NULL;
    }
    return cml_ones(M->row_n, M->col_n);
}


Matrix* cml_identity(int d) {
    if(d<=0){
        cml_error("Dimension d should be positive.");
        return NULL;
    }
    return basicIdentity(d);
}


Matrix* cml_transpose(Matrix* M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return basicTranspose(M);
}


Matrix* cml_add(Matrix* a, Matrix *b) {
    if(a==NULL || b==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(a->m==NULL || b->m==NULL){
        cml_error("Matrix-m is NULL.");
        return NULL;
    }
    if(a->row_n==b->row_n && a->col_n==b->col_n){
        return basicAdd(a, b);
    }else{
        cml_error("a and b have different shape.");
        return NULL;
    }
}


Matrix* cml_minus(Matrix* a, Matrix* b) {
    if(a==NULL || b==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(a->m==NULL || b->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return basicMinus(a, b);
}


Matrix* cml_numProd(double a, Matrix *M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return basicNumProd(a, M);
}


Matrix* cml_dot(Matrix *a, Matrix *b) {
    if(a==NULL || b==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(a->m==NULL || b->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(a->col_n!=b->row_n){
        cml_error("Cannot perform dot product on col(A) != row(B).");
        return NULL;
    }
    return basicDot(a, b);
}


double cml_det(Matrix* M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return basicDeterminant(M);
}


Matrix* cml_inv(Matrix *M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(basicDeterminant(M)==0){
        cml_error("Matrix is singular.");
        // warning or error ?
        return NULL;
    }
    return basicInverse(M);
}


Matrix* cml_diag(double *m, uint dim) {
    if(m==NULL){
        cml_error("Elements should be provided through m.");
        return NULL;
    }
    if(dim==0){
        cml_error("Dimension accepts only positive, dim=%d given.", dim);
        return NULL;
    }
    return basicDiag(m, dim);
}


Matrix* cml_cross_3d(Matrix *a, Matrix *b) {
    if(a==NULL || b==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(a->m==NULL || b->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }

    if(a->row_n*a->col_n==3 && b->row_n*b->col_n==3){
        return basic_3d_cross(a, b);
    }else{
        cml_error("A and B have to be 1x3 or 3x1.");
        return NULL;
    }
}


Matrix* cml_slice(Matrix *M, Matrix *R, Matrix *C) {
    if(M==NULL || R==NULL || C==NULL){
        cml_error("Matrix, Row_Idx or Col_Idx is NULL.");
        return NULL;
    }
    if(M->m==NULL || R->m==NULL || C->m==NULL){
        cml_error("Matrix, Row_Idx or Col_Idx->m is NULL.");
        return NULL;
    }
    if( (R->row_n!=1 && R->col_n!=1) ||
        (C->row_n!=1&&C->col_n!=1)
    ) {
        cml_error("R or C should be 1xn or nx1.");
        return NULL;
    }
    Matrix *r_tmp=NULL, *c_tmp=NULL;
    r_tmp = __cml_reformIndex__(M->row_n, R);
    c_tmp = __cml_reformIndex__(M->col_n, C);
    for(int i=0;i<r_tmp->col_n;i++){
        if(r_tmp->m[i]>=M->row_n){
            cml_error("row index %d out of range %d.", (int)r_tmp->m[i], M->row_n);
            return NULL;
        }
    }
    for(int i=0;i<c_tmp->col_n;i++){
        if(c_tmp->m[i]>=M->col_n){
            cml_error("col index %d out of range %d.", (int)c_tmp->m[i], M->col_n);
            return NULL;
        }
    }
    Matrix* Slice = basicSlice(M, r_tmp, c_tmp);
    basicFreeMatrix(&r_tmp);
    basicFreeMatrix(&c_tmp);
    return Slice;
}


Matrix* cml_rowSlice(Matrix *M, Matrix *R_Idx) {
    return cml_slice(M, R_Idx, basicRange(0,M->col_n,1));
}


Matrix* cml_colSlice(Matrix *M, Matrix *C_Idx) {
    return cml_slice(M, basicRange(0,M->row_n,1), C_Idx);
}


Matrix* cml_sum(Matrix* M, int axis) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(axis==-1 || axis==0 || axis==1){
        return basicSum(M, axis);
    }else{
        cml_error("axis only accepts -1, 0 or 1, %d received.", axis);
        return NULL;
    }
}


Matrix* cml_mean(Matrix *M, int axis) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(axis==-1 || axis==0 || axis==1){
        return basicMean(M, axis);
    }else{
        cml_error("axis only accepts -1, 0 or 1, %d received.", axis);
        return NULL;
    }
}


Matrix* cml_power(Matrix *M, int p) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(M->row_n!=M->col_n){
        cml_error("Only square matrix is available for power.");
        return NULL;
    }
    Matrix *re = NULL;
    if(p<0){
        if(basicDeterminant(M)==0){
            cml_error("Matrix is singular.");
            return NULL;
        }
        p = -p;
        if(p==p/2*2){
            // p is even
            re = basicCopy(M);
        }else{
            // p is odd
            re = basicInverse(M);
        }
    }else if(p==0){
        re = basicIdentity(M->row_n);
    }else{
        re = basicMatPow(M, p);
    }
    return re;
}


Matrix* cml_power_elementWise(Matrix *M, double p) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return basicElemPow(M, p);
}


Matrix* cml_sortVec(Matrix *vec) {
    if(vec==NULL){
        cml_error("Pointer vec is NULL.");
        return NULL;
    }
    if(vec->m==NULL){
        cml_error("vec->m is NULL.");
        return NULL;
    }
    return basicVecSort(vec);
}


Matrix* cml_sort(Matrix *M, int axis, int index) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(axis!=0 && axis!=1){
        cml_error("Axis accepts only 0 or 1.");
        return NULL;
    }
    if(index>=(axis==0?M->col_n:M->row_n)
      && index<-(axis==0?M->col_n:M->row_n)){
        cml_error("sort_m: index out of range.");
        return NULL;
    }
    if(index>=0){
        ;
    } else {
        index = (axis == 0 ? M->col_n : M->row_n) + index;
    }
    return basicMatSort(M, axis, index);
}


Matrix* cml_reverse(Matrix *vec) {
    if(vec==NULL){
        cml_error("Pointer vec is NULL.");
        return NULL;
    }
    if(vec->m==NULL){
        cml_error("vec->m is NULL.");
        return NULL;
    }
    return basicReverse(vec);
}


double cml_vecMax(Matrix *vec){
    if(vec==NULL){
        cml_error("Pointer vec is NULL.");
        return NULL;
    }
    if(vec->m==NULL){
        cml_error("vec->m is NULL.");
        return NULL;
    }
    return basicVecMax(vec);
}


double cml_vecMin(Matrix *vec){
    if(vec==NULL){
        cml_error("vec is NULL.");
        return NULl;
    }
    if(vec->m==NULL) {
        cml_error("vec->m is NULL.");
        return NULL;
    }
    return basicVecMin(vec);
}


Matrix* cml_max(Matrix *M, int axis) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(axis!=-1 && axis!=0 && axis !=1){
        cml_error("Axis should be -1, 0 or 1.");
        return NULL;
    }
    return basicMatMax(M, axis);
}


Matrix* cml_min(Matrix *M, int axis) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(axis!=-1 && axis!=0 && axis !=1){
        cml_error("Axis should be (-1, 0, 1).");
        return NULL;
    }
    return basicMatMin(M, axis);
}


int cml_apply(Matrix *From, Matrix *To, Matrix *Row_idx, Matrix*Col_idx) {
    if( !(From && To && Row_idx && Col_idx) ){
        cml_error("From, To, Row_idx, Col_idx cannot be NULL.");
        return ERROR;
    }
    if(
        !(From->m && To->m && Row_idx->m && Col_idx->m)
    ){
        cml_error("NULL found in Matrix->m(s).");
        return ERROR;
    }
    if(Row_idx->row_n!=1){
        Row_idx = basicTranspose(Row_idx);
    }
    if(Col_idx->row_n!=1){
        Col_idx = basicTranspose(Col_idx);
    }
    if(
        Row_idx->col_n > From->row_n ||
        Col_idx->col_n > From->col_n
    ){
        cml_error("Length of Row_idx or Col_idx is large than size of From.");
        return ERROR;
    }

    if(
        basicMatMax(Row_idx, -1)->m[0] > To->row_n ||
        basicMatMax(Col_idx, -1)->m[0] > To->col_n
    ) {
        cml_error("Index out of range. Value(s) of Row_idx or Col_idx is (are) larger than size of To.");
        return ERROR;
    }
    return basicApply(From, To, Row_idx, Col_idx);
}


Matrix* cml_shape(Matrix *M) {
    if(M==NULL) {
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_warning("Matrix->m is NULL.");
    }
    return basicShape(M);
}


Matrix* cml_reshape(Matrix *M, Matrix *Shape) {
    if(M==NULL || Shape==NULL){
        cml_error("Matrix or Shape is NULL.");
        return NULL;
    }
    if(Shape->m==NULL){
        cml_error("Shape->m is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_warning("Matrix->m is NULL.");
    }
    
    if(
        M->row_n * M->col_n != Shape->m[0] * Shape->m[1] ||
        Shape->m[0]<=0 || Shape->m[1]<=0
    ){
        cml_error("Incorrect Shape.");
        return NULL;
    }
    return basicReshape(M, Shape);
}


Matrix* cml_reshape2(Matrix *M, int row, int col) {
    Matrix *shape = basicZeros(1, 2);
    shape->m[0] = (double)row;
    shape->m[1] = (double)col;
    return cml_reshape(M, shape)
}


Matrix* cml_flatten(Matrix *M) {
    return = cml_reshape2(M, 1, M->row_n*M->col_n);
}


Matrix* cml_concatenate(Matrix *A, Matrix *B, int axis) {
    if (A==NULL || B==NULL) {
        cml_error("A or B is NULL.");
        return NULL;
    }
    if (A->m==NULL || B->m==NULL) {
        cml_error("A->m or B->m is NULL.");
        return NULL;
    }
    if (!(axis==0 || axis==1)) {
        cml_error("Axis can only be 0 for row wise or 1 for col wise concatenation.");
        return NULL;
    }
    return basicConcatenate(A, B, axis);
}

#endif
