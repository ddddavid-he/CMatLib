#ifndef MATRIX_CALCULATION
#define MATRIX_CALCULATION

#include <stdio.h>
#include <stdlib.h>

#include "basicCalculation.c"

#include "matrixCalculation.h"

#define ERROR 0
#define SUCCEED 1
#define WARN -1

#define TRUE 1
#define FALSE 0

typedef char bool;


int freeMatrix(Matrix **Mp) {
    // CAUTION: 
    // Mp have to be Pointer of Matrix Pointer.
    if(Mp==NULL){
        printf("freeMatrix: pointer of Matrix pointer Mp is NULL.\n");
        return ERROR;
    }else if((*Mp)==NULL){
        printf("freeMatrix: Matrix is already NULL.\n");
        return WARN;
    }else{
        basicFreeMatrix(Mp);
        return SUCCEED;
    }
}


int fillMatrixFromArray(double *array, int row, int col, Matrix *M) {
    if(array==NULL){
        printf("fillMatrixFromArray: array is NULL.\n");
        return ERROR;
    }else if(row<=0||col<=0){
        printf("fillMatrixFromArray: row and col should be positive.\n");
        return ERROR;
    }else if(M==NULL){
        printf("fillMatrixFromArray: Matrix is NULL.\n");
        return ERROR;
    }else{
        basicFillMatrixFromArray(array, row, col, M);
        return SUCCEED;
    }
}


Matrix* array2Matrix(double *array, int row, int col) {
    if(array==NULL){
        printf("array2Matrix: array is NULL.\n");
        return NULL;
    }else if(row<=0||col<=0){
        printf("array2Matrix: row and col should be positive.\n");
        return NULL;
    }else{
        return basicArray2Matrix(array, row, col);
    }
}


Matrix* __reformIndex__(int RowCol, Matrix *Indices) {
    Matrix *in = NULL;
    if(Indices->row_n!=1 && Indices->col_n==1){
        in = basicTranspose(Indices);
    }else{
        in = basicCopy(Indices);
    }

    for(int i=0;i<in->row_n;i++){
        for(int j=0;j<Indices->col_n;j++){
            int index = giveIndex(in, i, j);
            if(in->m[index]<0){
                if(in->m[index]<-RowCol){
                    printf("__reformIndex__: index %d out of range %d.\n", (int)in->m[index], RowCol);
                    exit(-1);
                }
                in->m[index] += RowCol;
            }
        }
    }
    return in;
}


int getIndex(Matrix* M, int i, int j) {
    if(M==NULL){
        printf("getIndex: Matrix is NULL.\n");
        exit(-1);
    }else{
        if(M->m==NULL){
            printf("W: getIndex: Matrix->m is NULL.\n");
        }
    }
    int row=M->row_n, col=M->col_n;
    if(i>=row && i<-row){
        printf("getIndex: row index i out of range. Row=%d\n", row);
        exit(-1);
    }else if(j>=col && j<-col){
        printf("getIndex: column index j out of range. Col=%d\n", col);
        exit(-1);
    }else{
        if(i<0){
            i += row;
        }
        if(j<0){
            j += col;
        }
    }
    return giveIndex(M, i, j);
}


int showMatrix(Matrix *M) {
    if(M==NULL){
        printf("showMatrix: Matrix is NULL\n");
        return ERROR;
    }else{
        if(M->m==NULL){
            printf("showMatrix: Matrix->m is NULL\n");
            return ERROR;
        }
    }
    basicShowMatrix(M);
    return SUCCEED;
}


Matrix* copy_m(Matrix *M) {
    Matrix *re = malloc(sizeof(Matrix));
    if(re==NULL){
        printf("copy_m: Fail to malloc.\n");
        return NULL;
    }else if(M==NULL){
        printf("copy_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("copy_m: Matrix->m is NULL.\n");
        return NULL;
    }else{
        return basicCopy(M);
    }
}


Matrix* range(int from, int to, int skip) {
    if(from==to){
        if(skip==0){
            return basicNumProd(from, basicOnes(1, 1));
        }else{
            printf("range: Empty range \n");
            return NULL;
        }
    }else if(from<to){
        if(skip==0){
            skip = 1;
        }else if(skip<0){
            printf("range: Emtpy range.\n");
            return NULL;
        }
        return basicRange(from, to, skip);
    }else{ 
        // from > to
        if(skip==0){
            skip = -1;
        }else if(skip>0){
            printf("range: Empty range.\n");
            return NULL;
        }
        return basicRange(from, to, skip);
    }
}


Matrix* zeros_m(int row, int col) {
    if(row*col<=0){
        printf("zeros_m: Cannot create a Matrix with non-positive row or column number.\n");
        return NULL;
    }else{
        return basicZeros(row, col);
    }
}


Matrix* ones_m(int row, int col) {
    if(row*col<=0){
        printf("ones_m: Cannot create a Matrix with non-positive row or column number.\n");
        return NULL;
    }else{
        return basicOnes(row, col);
    }
}


Matrix* zerosLike_m(Matrix* M) {
    if(M==NULL){
        printf("zerosLike_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("W: zerosLike_m: Matrix->m is NULL\n");
        return zeros_m(M->row_n, M->col_n);
    }else{
        return zeros_m(M->row_n, M->col_n);
    } 
}


Matrix* onesLike_m(Matrix* M) {
    if(M==NULL){
        printf("onesLike_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("W: onesLike_m: Matrix->m is NULL\n");
        return ones_m(M->row_n, M->col_n);
    }else{
        return ones_m(M->row_n, M->col_n);
    }
}


Matrix* identity_m(int d) {
    if(d<=0){
        printf("identity_m: Dimension d should be positive.\n");
        return NULL;
    }else{
        return basicIdentity(d);
    }
}


Matrix* identityLike_m(Matrix *M) {
    if(M==NULL){
        printf("identityLike_m: Matrix is NULL\n");
        return NULL;
    }else if(M->m==NULL){
        printf("W: identityLike_m: Matrix->m is NULL\n");
        return NULL;
    }else{
        if(M->row_n!=M->col_n){
            printf("identityLike_m: Only square matrix is acceptable for this function.\n");
            return NULL;
        }else{
            return basicIdentityLike(M);
        }
    }
}


Matrix* trans_m(Matrix* M) {
    if(M==NULL){
        printf("trans_m: Matrix is NULL\n");
        return NULL;
    }else{
        if(M->m==NULL){
            printf("trans_m: Matrix->m is NULL\n");
            return NULL;
        }else{
            return basicTranspose(M);
        }
    }
}


Matrix* add_m(Matrix* a, Matrix *b) {
    if(a==NULL || b==NULL){
        printf("add_m: Matrix is NULL.\n");
        return NULL;
    }else{
        if(a->m==NULL || b->m==NULL){
            printf("add_m: Matrix-m is NULL.\n");
            return NULL;
        }
        if(a->row_n==b->row_n && a->col_n==b->col_n){
            return basicAdd(a, b);
        }else{
            printf("add_m: a and b have different shape.\n");
            return NULL;
        }
    }
}


Matrix* nProd_m(double a, Matrix *M) {
    if(M==NULL){
        printf("nProd_m: Matrix is NULL.\n");
        return NULL;
    }else{
        if(M->m==NULL){
            printf("nProd_m: Matrix->m is NULL.\n");
            return NULL;
        }else{
            return basicNumProd(a, M);
        }
    }
}


Matrix* minus_m(Matrix* a, Matrix* b) {
    if(a==NULL || b==NULL){
        printf("minus_m: Matrix is NULL.\n");
        return NULL;
    }else{
        if(a->m==NULL || b->m==NULL){
            printf("minus_m: Matrix->m is NULL.\n");
            return NULL;
        }else{
            return basicMinus(a, b);
        }
    }
}


Matrix* dot_m(Matrix *a, Matrix *b) {
    if(a==NULL || b==NULL){
        printf("dot_m: Matrix is NULL.\n");
        return NULL;
    }else if(a->m==NULL || b->m==NULL){
        printf("dot_m: Matrix->m is NULL.\n");
        return NULL;
    }else if(a->col_n!=b->row_n){
        printf("dot_m: Cannot calculate dot product of col(A) != row(B).\n");
        return NULL;
    }else{
        return basicDot(a, b);
    }
}


double det_m(Matrix* M) {
    if(M==NULL){
        printf("det_m: Matrix is NULL.\n");
        exit(-1);
    }else if(M->m==NULL){
        printf("det_m: Matrix->m is NULL.\n");
        exit(-1);
    }else{
        return basicDeterminant(M);
    }
}


Matrix* inv_m(Matrix *M) {
    if(M==NULL){
        printf("inv_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("inv_m: Matrix->m is NULL.\n");
        return NULL;
    }else{
        if(basicDeterminant(M)==0){
            printf("det_m: Matrix is singular.\n");
            return NULL;
        }else{
            return basicInverse(M);
        }
    }
}


Matrix* diag_m(double *m, int dim) {
    if(m==NULL){
        printf("diag_m: elem is NULL.\n");
        return NULL;
    }else if(dim<=0){
        printf("diag_m: Dimension accepts only positive, dim=%d given.\n", dim);
        return NULL;
    }else{
        return basicDiag(m, dim);
    }
}


Matrix* cross_3d_m(Matrix *a, Matrix *b) {
    if(a==NULL || b==NULL){
        printf("cross_3d_m: Matrix is NULL.\n");
        return NULL;
    }else if(a->m==NULL || b->m==NULL){
        printf("cross_3d_m: Matrix->m is NULL.\n");
        return NULL;
    }

    if(a->row_n*a->col_n==3 && b->row_n*b->col_n==3){
        return basic_3d_cross(a, b);
    }else{
        printf("cross_3d_m: a and b have to be 1x3 or 3x1.\n");
        return NULL;
    }
}


Matrix* slice_m(Matrix *M, Matrix *R, Matrix *C) {
    if(M==NULL || R==NULL || C==NULL){
        printf("slice_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL || R->m==NULL || C->m==NULL){
        printf("slice_m: Matrix->m is NULL.\n");
        return NULL;
    }else if((R->row_n!=1&&R->col_n!=1)
          || (C->row_n!=1&&C->col_n!=1)){
        printf("slice_m: R or C should be 1xn or nx1.\n");
        return NULL;
    }else{
        Matrix *r_tmp=NULL, *c_tmp=NULL;
        r_tmp = __reformIndex__(M->row_n, R);
        c_tmp = __reformIndex__(M->col_n, C);
        for(int i=0;i<r_tmp->col_n;i++){
            if(r_tmp->m[i]>=M->row_n){
                printf("slice_m: row indices %d out of range %d.\n", (int)r_tmp->m[i], M->row_n);
                return NULL;
            }
        }
        for(int i=0;i<c_tmp->col_n;i++){
            if(c_tmp->m[i]>=M->col_n){
                printf("slice_m: col indices %d out of range %d.\n", (int)c_tmp->m[i], M->col_n);
                return NULL;
            }
        }
        Matrix* Slice = basicSlice(M, r_tmp, c_tmp);
        basicFreeMatrix(&r_tmp);
        basicFreeMatrix(&c_tmp);
        return Slice;
    }
}


Matrix* rowSlice_m(Matrix *M, Matrix *RS) {
    Matrix *re = slice_m(M, RS, basicRange(0,M->col_n,1));
    if(!re) {
        printf("rowSlice_m: Error. See message of slice_m above.\n");
    }
    return re;
}


Matrix* colSlice_m(Matrix *M, Matrix *CS) {
    Matrix *re = slice_m(M, basicRange(0,M->row_n,1), CS);
    if(!re) {
        printf("colSlice_m: Error. See message of slice_m above.\n");
    }
    return re;
}


Matrix* sum_m(Matrix* M, int axis) {
    if(M==NULL){
        printf("sum_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("sum_m: Matrix->m is NULL.\n");
        return NULL;
    }
    if(axis==-1 || axis==0 || axis==1){
        return basicSum(M, axis);
    }else{
        printf("sum_m: axis not avaliable.\n");
        return NULL;
    }
}


Matrix* mean_m(Matrix *M, int axis) {
    if(M==NULL){
        printf("mean_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("mean_m: Matrix->m is NULL.\n");
        return NULL;
    }
    if(axis==-1 || axis==0 || axis==1){
        return basicMean(M, axis);
    }else{
        printf("mean_m: axis not avaliable.\n");
        return NULL;
    }
}


Matrix* matPow_m(Matrix *M, int p) {
    if(M==NULL){
        printf("matPow_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("matPow_m: Matrix->m is NULL.\n");
        return NULL;
    }else if(M->row_n!=M->col_n){
        printf("matPow_m: Matrix is not a square matrix.\n");
        return NULL;
    }
    Matrix *re = NULL;
    if(p<0){
        if(basicDeterminant(M)==0){
            printf("matPow_m: Matrix is singular.\n");
            return NULL;
        }
        p = -p;
        if(p==p/2*2){
            re = basicCopy(M);
        }else{
            re = basicInverse(M);
        }
    }else if(p==0){
        re = basicIdentityLike(M);
    }else{
        re = basicMatPow(M, p);
    }
    return re;
}


Matrix* elemPow_m(Matrix* M, double p) {
    if(M==NULL){
        printf("elemPow_m: Matrix is NULL.\n");
        return NULL;
    }else if(M->m==NULL){
        printf("elemPow_m: Matrix->m is NULL.\n");
        return NULL;
    }else{
        return basicElemPow(M, p);
    }
}


Matrix* sortVec_m(Matrix *vec) {
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


Matrix* sort_m(Matrix *M, int axis, int index) {
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


Matrix* reverse_m(Matrix *vec) {
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


double vecMax_m(Matrix *vec){
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


double vecMin_m(Matrix *vec){
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


Matrix* max_m(Matrix *M, int axis) {
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


Matrix* min_m(Matrix *M, int axis) {
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


int apply_m(Matrix *From, Matrix *To, Matrix *RP, Matrix*CP) {
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


Matrix* shape_m(Matrix *M) {
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


Matrix* reshape_m(Matrix *M, Matrix *Shape) {
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


Matrix* reshape2_m(Matrix *M, int row, int col) {
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


Matrix* flatten_m(Matrix *M) {
    Matrix *re = reshape2_m(M, 1, M->row_n*M->col_n);
    if(!re){
        printf("flatten_m: Error. See reshape_m above.\n");
        return NULL;
    }
    return re;
}


#endif
