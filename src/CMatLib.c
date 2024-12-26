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
    char* type_name = (char*) calloc(7, sizeof(char));
    if(type==0){
        type_name = "Error";
    } else {
        type_name = "Warning";
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





int cml_freeMatrix(cml_Matrix_t **Mp) {
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

    cml_basicFreeMatrix(Mp);
    return SUCCEED;
}


cml_Matrix_t* cml_empty(uint row, uint col) {
    cml_Matrix_t* M = (cml_Matrix_t*) malloc(sizeof(cml_Matrix_t));
    M->m = calloc(row*col, sizeof(double));
    return M;
}


int cml_fillMatrixFromArray(double *array, cml_Matrix_t *M) {
    // fill values from array to matrix M
    // notice: Matrix->m should be allocated before this function
    if( array==NULL ) {
        cml_error("array is NULL.");
        return ERROR;
    }
    if( M->row_n<=0 || M->col_n<=0 ) {
        cml_error("row and col should be positive.");
        return ERROR;
    }
    if( M==NULL ) {
        cml_error("Matrix is NULL.");
        return ERROR;
    }
    if( M->m == NULL ) {
        cml_error("Matrix->m is not allocated yet.");
        return ERROR;
    }

    cml_basicFillMatrixFromArray(array, M);
    return SUCCEED;
}


cml_Matrix_t* cml_arrayToMatrix(double *array, int row, int col) {
    // create a Matrix from values of array
    cml_Matrix_t *M = (cml_Matrix_t*) malloc(sizeof(cml_Matrix_t));
    M-> m = (double*) calloc(row*col, sizeof(double));
    cml_fillMatrixFromArray(array, M);
    return M;
}


cml_Matrix_t* __cml_reformIndex__(int RowCol, cml_Matrix_t *Indices) {
    // enable index to be in (-l, l-1)
    // reform all the elements in Indices by RowCol (l).
    cml_Matrix_t *idx = cml_basicFlatten(Indices);
    cml_Matrix_t *shape = cml_basicShape(Indices);


    for(int i=0;i<idx->col_n;i++){
        if(idx->m[i] < -RowCol){
            cml_error("index %d out of range %d.", (int)idx->m[index], RowCol);
            return NULL;
        } else if (-RowCol <= idx->m[i] && idx->m[i] < 0 ) {
            idx->m[i] += RowCol;
        } else if (idx->m[i] >= 0 && idx->m[i] < RowCol) {
            ;
        } else {
            cml_error("index %d out of range %d.", (int)idx->m[index], RowCol);
            return NULL;
        }
    }
    cml_Matrix_t *tmp = cml_basicReshape(idx, shape);
    cml_basicFreeMatrix(&idx);
    idx = tmp;
    return idx;
}


int cml_getIndex(cml_Matrix_t* M, int row, int col) {
    // transform row index and col index into 1D index
    // usage:
    // M->m[cml_getIndex(M, 2, 2)]
    if(M==NULL){
        cml_error("Matrix is NULL.");
    }
    if(M->m==NULL){
        cml_warning("Matrix->m is NULL.");
    }
    if(row >= M->row_n && row < -M->row_n){
        cml_error("getIndex: row index i out of range. Row=%d", row);
        return ERROR;
    }else if(col >= M->col_n && col < -M->col_n){
        cml_error("getIndex: column index j out of range. Col=%d", col);
        return ERROR;
    }

    if(row<0){
        row += M->row_n;
    }
    if(col<0){
        col += M->col_n;
    }
    return cml_basicGetIndex(M, row, col);
}


int cml_showMatrix(cml_Matrix_t *M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return ERROR;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL");
        return ERROR;
    }
    cml_basicShowMatrix(M);
    return SUCCEED;
}


cml_Matrix_t* cml_copy(cml_Matrix_t *M) {
    if( M==NULL || M->m==NULL ){
        cml_error("Matrix is NULL.");
        return NULL;
    }

    return cml_basicCopy(M);
}


cml_Matrix_t* cml_range(int from, int to, int step) {
    if(from==to){
        if(step==0){
            return cml_basicNumProd(from, cml_basicOnes(1, 1));
        }else{
            cml_error("Empty range.");
            return NULL;
        }
    }else if(from<to){
        if(step<=0){
            cml_error("Emtpy range.");
            return NULL;
        }
        return cml_basicRange(from, to, step);
    }else{ 
        // from > to
        if(step>=0){
            cml_error("Empty range.");
            return NULL;
        }
        return cml_basicRange(from, to, step);
    }
}


cml_Matrix_t* cml_zeros(int row, int col) {
    if(row*col<=0){
        cml_error("Cannot create a Matrix with non-positive row or column number.");
        return NULL;
    }
    return cml_basicZeros(row, col);
}


cml_Matrix_t* cml_ones(int row, int col) {
    if(row*col<=0){
        cml_error("Cannot create a Matrix with non-positive row or column number.");
        return NULL;
    }
    return cml_basicOnes(row, col);
}


cml_Matrix_t* cml_zerosLike(cml_Matrix_t* M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    return cml_zeros(M->row_n, M->col_n);
}


cml_Matrix_t* cml_onesLike(cml_Matrix_t* M) {
    if(M==NULL) {
        cml_error("Matrix is NULL.");
        return NULL;
    }
    return cml_ones(M->row_n, M->col_n);
}


cml_Matrix_t* cml_identity(int d) {
    if(d<=0){
        cml_error("Dimension d should be positive.");
        return NULL;
    }
    return cml_basicIdentity(d);
}


cml_Matrix_t* cml_transpose(cml_Matrix_t* M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return cml_basicTranspose(M);
}


cml_Matrix_t* cml_add(cml_Matrix_t* a, cml_Matrix_t *b) {
    if(a==NULL || b==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(a->m==NULL || b->m==NULL){
        cml_error("Matrix-m is NULL.");
        return NULL;
    }
    if(a->row_n==b->row_n && a->col_n==b->col_n){
        return cml_basicAdd(a, b);
    }else{
        cml_error("a and b have different shape.");
        return NULL;
    }
}


cml_Matrix_t* cml_minus(cml_Matrix_t* A, cml_Matrix_t* B) {
    if(A==NULL || B==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(A->m==NULL || B->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return cml_basicMinus(A, B);
}


cml_Matrix_t* cml_numProd(double a, cml_Matrix_t *M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return cml_basicNumProd(a, M);
}


cml_Matrix_t* cml_dot(cml_Matrix_t *a, cml_Matrix_t *b) {
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
    return cml_basicDot(a, b);
}


double cml_det(cml_Matrix_t* M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
//        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
//        return NULL;
    }
    return cml_basicDeterminant(M);
}


cml_Matrix_t* cml_inv(cml_Matrix_t *M) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(cml_basicDeterminant(M)==0){
        cml_error("Matrix is singular.");
        // warning or error ?
        return NULL;
    }
    return cml_basicInverse(M);
}


cml_Matrix_t* cml_diag(double *m, int dim) {
    if(m==NULL){
        cml_error("Elements should be provided through m.");
        return NULL;
    }
    if(dim==0){
        cml_error("Dimension accepts only positive, dim=%d given.", dim);
        return NULL;
    }
    return cml_basicDiag(m, dim);
}


cml_Matrix_t* cml_cross(cml_Matrix_t *a, cml_Matrix_t *b) {
    if(a==NULL || b==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(a->m==NULL || b->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }

    if(a->row_n*a->col_n==3 && b->row_n*b->col_n==3){
        return cml_basicCross(a, b);
    }else{
        cml_error("A and B have to be 1x3 or 3x1.");
        return NULL;
    }
}


cml_Matrix_t* cml_slice(cml_Matrix_t *M, cml_Matrix_t *R, cml_Matrix_t *C) {
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
    cml_Matrix_t *r_tmp=NULL, *c_tmp=NULL;
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
    cml_Matrix_t* Slice = cml_basicSlice(M, r_tmp, c_tmp);
    cml_basicFreeMatrix(&r_tmp);
    cml_basicFreeMatrix(&c_tmp);
    return Slice;
}


cml_Matrix_t* cml_rowSlice(cml_Matrix_t *M, cml_Matrix_t *R_Idx) {
    return cml_slice(M, R_Idx, cml_basicRange(0,M->col_n,1));
}


cml_Matrix_t* cml_colSlice(cml_Matrix_t *M, cml_Matrix_t *C_Idx) {
    return cml_slice(M, cml_basicRange(0,M->row_n,1), C_Idx);
}


cml_Matrix_t* cml_sum(cml_Matrix_t* M, int axis) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(axis==-1 || axis==0 || axis==1){
        return cml_basicSum(M, axis);
    }else{
        cml_error("axis only accepts -1, 0 or 1, %d received.", axis);
        return NULL;
    }
}


cml_Matrix_t* cml_mean(cml_Matrix_t *M, int axis) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    if(axis==-1 || axis==0 || axis==1){
        return cml_basicMean(M, axis);
    }else{
        cml_error("axis only accepts -1, 0 or 1, %d received.", axis);
        return NULL;
    }
}


cml_Matrix_t* cml_power(cml_Matrix_t *M, int p) {
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
    cml_Matrix_t *re = NULL;
    if(p<0){
        if(cml_basicDeterminant(M)==0){
            cml_error("Matrix is singular.");
            return NULL;
        }
        p = -p;
        if(p==p/2*2){
            // p is even
            re = cml_basicCopy(M);
        }else{
            // p is odd
            re = cml_basicInverse(M);
        }
    }else if(p==0){
        re = cml_basicIdentity(M->row_n);
    }else{
        re = cml_basicMatPow(M, p);
    }
    return re;
}


cml_Matrix_t* cml_power_elementWise(cml_Matrix_t *M, double p) {
    if(M==NULL){
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_error("Matrix->m is NULL.");
        return NULL;
    }
    return cml_basicElemPow(M, p);
}


cml_Matrix_t* cml_sortVec(cml_Matrix_t *vec) {
    if(vec==NULL){
        cml_error("Pointer vec is NULL.");
        return NULL;
    }
    if(vec->m==NULL){
        cml_error("vec->m is NULL.");
        return NULL;
    }
    return cml_basicVecSort(vec);
}


cml_Matrix_t* cml_sort(cml_Matrix_t *M, int axis, int index) {
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
    return cml_basicMatSort(M, axis, index);
}


cml_Matrix_t* cml_reverse(cml_Matrix_t *vec) {
    if(vec==NULL){
        cml_error("Pointer vec is NULL.");
        return NULL;
    }
    if(vec->m==NULL){
        cml_error("vec->m is NULL.");
        return NULL;
    }
    return cml_basicReverse(vec);
}


double cml_vecMax(cml_Matrix_t *vec){
    if(vec==NULL){
        cml_error("Pointer vec is NULL.");
//        return NULL;
    }
    if(vec->m==NULL){
        cml_error("vec->m is NULL.");
//        return NULL;
    }
    return cml_basicVecMax(vec);
}


double cml_vecMin(cml_Matrix_t *vec){
    if(vec==NULL){
        cml_error("vec is NULL.");
//        return NULl;
    }
    if(vec->m==NULL) {
        cml_error("vec->m is NULL.");
//        return NULL;
    }
    return cml_basicVecMin(vec);
}


cml_Matrix_t* cml_max(cml_Matrix_t *M, int axis) {
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
    return cml_basicMatMax(M, axis);
}


cml_Matrix_t* cml_min(cml_Matrix_t *M, int axis) {
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
    return cml_basicMatMin(M, axis);
}


int cml_apply(cml_Matrix_t *From, cml_Matrix_t *To, cml_Matrix_t *Row_idx, cml_Matrix_t*Col_idx) {
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
    int length_of_RI = (int) Row_idx->row_n * Row_idx->col_n;
    int length_of_CI = (int) Col_idx->row_n * Col_idx->col_n;
    int length_of_From = (int) From->row_n * From->col_n;
    if( length_of_RI != length_of_CI ){
        cml_error("Length of Row_idx and Col_idx does not match. Same length is required.");
    }
    if( length_of_RI > length_of_From ){
        cml_error("Length of Row_idx and Col_idx is large than size of From.");
        return ERROR;
    }
    if(
        cml_basicMatMax(Row_idx, -1)->m[0] > To->row_n ||
        cml_basicMatMax(Col_idx, -1)->m[0] > To->col_n
    ) {
        cml_error("Index out of range. Value(s) of Row_idx or Col_idx is (are) larger than size of To.");
        return ERROR;
    }
    return cml_basicApply(From, To, Row_idx, Col_idx);
}


cml_Matrix_t* cml_shape(cml_Matrix_t *M) {
    if(M==NULL) {
        cml_error("Matrix is NULL.");
        return NULL;
    }
    if(M->m==NULL){
        cml_warning("Matrix->m is NULL.");
    }
    return cml_basicShape(M);
}


cml_Matrix_t* cml_reshape(cml_Matrix_t *M, cml_Matrix_t *Shape) {
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
    return cml_basicReshape(M, Shape);
}


cml_Matrix_t* cml_reshape2(cml_Matrix_t *M, int row, int col) {
    cml_Matrix_t *shape = cml_basicZeros(1, 2);
    shape->m[0] = (double)row;
    shape->m[1] = (double)col;
    return cml_reshape(M, shape);
}


cml_Matrix_t* cml_flatten(cml_Matrix_t *M) {
    return cml_reshape2(M, 1, M->row_n*M->col_n);
}


cml_Matrix_t* cml_concatenate(cml_Matrix_t *A, cml_Matrix_t *B, int axis) {
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
    return cml_basicConcatenate(A, B, axis);
}

#endif
