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
        const char *type_name = "Error"
    } else {
        const char *type_name = "Warning"
    }
    fprintf(stderr, "%s in func %s: (File: %s, Line: %d) ", type_name, func, file, line);
    vfprintf(stderr, message, args);
    fprintf(stderr, "\n");
    va_end(args);
    if(type==0){
        exit(-1)
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
        cml_error("Matrix->m is not allocated yet.")
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
    return zeros_m(M->row_n, M->col_n);
}


Matrix* cml_onesLike(Matrix* M) {
    if(M==NULL) {
        cml_error("Matrix is NULL.");
        return NULL;
    }
    return ones_m(M->row_n, M->col_n);
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
    r_tmp = __reformIndex__(M->row_n, R);
    c_tmp = __reformIndex__(M->col_n, C);
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
    return slice_m(M, basicRange(0,M->row_n,1), C_Idx);
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
        cml_error("only square matrix is available for power.");
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


Matrix* cml_sortVec_m(Matrix *vec) {
    if(vec==NULL){
        printf("sortVec: Vec is NULL.\n");
        return NULL;
    }else if(vec->m==NULL){
        printf("sort_m: Vec->m is NULL.\n");
        return NULL;
    }else{
        return basicVecSort(vec);
    }
}


Matrix* cml_sort_m(Matrix *M, int axis, int index) {
    if(M==NULL){
        printf("sort_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("sort_m: Matrix->m is NULL.\n");
        return NULL;
    }else{
        if(axis!=0 && axis!=1){
            printf("sort_m: axis accepts only 0 or 1.\n");
            return NULL;
        }else if(index>=(axis==0?M->col_n:M->row_n) 
              && index<-(axis==0?M->col_n:M->row_n)){
            printf("sort_m: index out of range.\n");
            return NULL;
        }else{
            index = index>=0?index:
                (axis==0?M->col_n:M->row_n)+index;
            return basicMatSort(M, axis, index);
        }
    }
}


Matrix* cml_reverse_m(Matrix *vec) {
    if(vec==NULL){
        printf("reverse_m: Vec is NULL.\n");
        return NULL;
    }else if(vec->m==NULL){
        printf("reverse_m: Vec->m is NULL.\n");
        return NULL;
    }else{
        return basicReverse(vec);
    }
}


double cml_vecMax_m(Matrix *vec){
    if(vec==NULL){
        printf("vecMax_m: Vec is NULL.\n");
        exit(-1);
    }else if(vec->m==NULL){
        printf("vecMax_m: Vec->m is NULL.\n");
        exit(-1);
    }else{
        return basicVecMax(vec);
    }
}


double cml_vecMin_m(Matrix *vec){
    if(vec==NULL){
        printf("vecMin_m: Vec is NULL.\n");
        exit(-1);
    }else if(vec->m==NULL){
        printf("vecMin_m: Vec->m is NULL.\n");
        exit(-1);
    }else{
        return basicVecMin(vec);
    }
}


Matrix* cml_max_m(Matrix *M, int axis) {
    if(M==NULL){
        printf("max_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("max_m: Matrix->m is NULL.\n");
        return NULL;
    }else{
        if(axis!=-1 && axis!=0 && axis !=1){
            printf("max_m: axis should be (-1, 0, 1).\n");
            return NULL;
        }else{
            return basicMatMax(M, axis);
        }
    }
}


Matrix* cml_min_m(Matrix *M, int axis) {
    if(M==NULL){
        printf("min_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("min_m: Matrix->m is NULL.\n");
        return NULL;
    }else{
        if(axis!=-1 && axis!=0 && axis !=1){
            printf("min_m: axis should be (-1, 0, 1).\n");
            return NULL;
        }else{
            return basicMatMin(M, axis);
        }
    }
}


int cml_apply_m(Matrix *From, Matrix *To, Matrix *RP, Matrix*CP) {
    if(!(From&&To&&RP&&CP)){
        printf("apply_m: At least one of the Matrixes is NULL.\n");
        return ERROR;
    }else if(!(From->m&&To->m&&RP->m&&CP->m)){
        printf("apply_m: At least one of Matrix->m(s) is NULL.\n");
        return ERROR;
    }else{
        if(RP->row_n!=1){
            RP = basicTranspose(RP);
        }
        if(CP->row_n!=1){
            CP = basicTranspose(CP);
        }
        if(RP->col_n>From->row_n||CP->col_n>From->col_n){
            printf("apply_m: Index out of range. Length of RP or CP is large than size of From.\n");
            return ERROR;
        }else if(RP->col_n<From->row_n||CP->col_n<From->col_n){
            printf("W: apply_m: Length of RP or CP does not match (smaller) size of From.\n");
            return WARN;
        }else{
            if(RP->col_n>To->row_n||CP->col_n>To->col_n){
                printf("apply_m: Index out of range. Length of RP or CP is larger than size of To.\n");
                return ERROR;
            }else{
                if(basicApply(From, To, RP, CP)){
                    return SUCCEED;
                }
                else{
                    return ERROR;
                }
            }
        }
    }
}


Matrix* cml_shape_m(Matrix *M) {
    if(M==NULL) {
        printf("shape_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("W: shape_m: Matrix->m is NULL.\n");
        return basicShape(M);
    }else{
        return basicShape(M);
    }
}


Matrix* cml_reshape_m(Matrix *M, Matrix *Shape) {
    if(M==NULL||Shape==NULL){
        printf("reshape_m: Matrix or Shape is NULL.\n");
        return NULL;
    }else if(Shape->m==NULL){
        printf("reshape_m: Shape->m is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("W: reshape_m: Matrix->m is NULL.\n");
    }
    
    if(M->row_n*M->col_n!=Shape->m[0]*Shape->m[1]
    || Shape->m[0]<=0 || Shape->m[1]<=0){
        printf("reshape_m: Shape and Matrix do not match.\n");
        return NULL;
    }else{
        return basicReshape(M, Shape);
    }
}


Matrix* cml_reshape2_m(Matrix *M, int row, int col) {
    Matrix *shape = basicZeros(1, 2);
    shape->m[0] = (double)row;
    shape->m[1] = (double)col;
    Matrix *re = reshape_m(M, shape);
    if(!re){
        printf("reshape2_m: Error. See reshape_m above.\n");
        return NULL;
    }
    basicFreeMatrix(&shape);
    return re;
}


Matrix* cml_flatten_m(Matrix *M) {
    Matrix *re = reshape2_m(M, 1, M->row_n*M->col_n);
    if(!re){
        printf("flatten_m: Error. See reshape_m above.\n");
        return NULL;
    }
    return re;
}


#endif
