#ifndef CML_BASIC
#define CML_BASIC

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basic.h" 



void cml_basicFreeMatrix(cml_Matrix_t **Mp) {
    free((*Mp)->m);
    free(*Mp);
    *Mp = NULL;
}


int cml_basicGetIndex(cml_Matrix_t *M, int i, int j) {
    return M->col_n * i + j;
}


void cml_basicShowMatrix(cml_Matrix_t *M) {
    int i,j,k;
    char format[11] = "-8.3lf  \0";

    format[1] = (char) 48 + NUM_LEN;
    format[3] = (char) 48 + FLOAT_LEN;
    // printf("%s\n", format);
    printf("\t_");
    for(k = 0;
        k < M->col_n*NUM_LEN + (M->col_n)*2 + 1;
        k++){
            printf(" ");
    }
    printf("_\n");
        
    for(i=0;i<M->row_n;i++){
        printf("\t|  ");
        for(j=0;j<M->col_n;j++){
            int index = cml_basicGetIndex(M, i, j);
            printf("%-8.3lf  ", M->m[index]);
        }
        printf("|\n");
    }

    printf("\t_");
    for(k = 0;
        k < M->col_n*NUM_LEN + (M->col_n)*2 + 1;
        k++){
            printf(" ");
    }
    printf("_\n");
}


void cml_basicFillMatrixFromArray(double* array, cml_Matrix_t* M) {
    int i, j, idx;
    const int row = M->row_n;
    const int col = M->col_n;
    double *m = M->m;
    for(i=0;i<row;i++) {
        for(j=0;j<col;j++) {
            idx = cml_basicGetIndex(M, i, j);
            m[idx] = array[idx];
        }
    }
}


cml_Matrix_t* cml_basicArray2Matrix(double* array, int row, int col) {
    cml_Matrix_t *M = malloc(sizeof(cml_Matrix_t));
    M->m = (double*) calloc(row * col, sizeof(double));
    cml_basicFillMatrixFromArray(array, M);
    return M;
}


cml_Matrix_t* cml_basicCopy(cml_Matrix_t *M) {
    cml_Matrix_t *result = malloc(sizeof(cml_Matrix_t));
    result->m = (double*) calloc(M->row_n*M->col_n, sizeof(double));
    result->row_n = M->row_n;
    result->col_n = M->col_n;
    for(int i=0;i<M->row_n*M->col_n;i++){
        result->m[i] = M->m[i];
    }
    return result;
}


cml_Matrix_t* cml_basicRange(int from, int to, int step) {
    const int length = (int)ceil(((double)(to-from))/((double)step));
    cml_Matrix_t * vec = cml_basicZeros(1, length);
    if(step>0){
        for(int i=from,j=0; i<to; i+=step,j++){
            vec->m[j] = i;
        }
    }else{
        for(int i=from,j=0; i>to; i+=step,j++){
            vec->m[j] = i;
        }
    }

    return vec;
} 


cml_Matrix_t* cml_basicEmpty(int row, int col) {
    cml_Matrix_t *M = (cml_Matrix_t*) malloc(sizeof(cml_Matrix_t));
    M->row_n = row;
    M->col_n = col;
    M->m = calloc(row*col, sizeof(double));
    return M;
}


cml_Matrix_t* cml_basicZeros(int row, int col) {
    cml_Matrix_t *M = malloc(sizeof(cml_Matrix_t));
    M->row_n = row;
    M->col_n = col;
    M->m = calloc(row*col, sizeof(double));
    for(int i=0;i<row*col;i++) {
        M->m[i] = 0.;
    }
    return M;
}


cml_Matrix_t* cml_basicOnes(int row, int col) {
    cml_Matrix_t * M = malloc(sizeof(cml_Matrix_t));
    M->row_n = row;
    M->col_n = col;
    M->m = calloc(row*col, sizeof(double));
    for(int i=0;i<row*col;i++) {
        M->m[i] = 1.;
    }
    return M;
}


cml_Matrix_t* cml_basicZerosLike(cml_Matrix_t *M) {
    return cml_basicZeros(M->row_n, M->col_n);
}


cml_Matrix_t* cml_basicOnesLike(cml_Matrix_t *M) {
    return cml_basicOnes(M->row_n, M->col_n);
}


cml_Matrix_t* cml_basicIdentity(int d) {
    cml_Matrix_t *re = cml_basicZeros(d, d);
    for(int i=0;i<d;i++){
        re->m[cml_basicGetIndex(re, i, i)] = 1.;
    }
    return re;
} 


cml_Matrix_t* cml_basicIdentityLike(cml_Matrix_t* M) {
    cml_Matrix_t *re = cml_basicIdentity(M->row_n);
    return re;
}


cml_Matrix_t* cml_basicTranspose(cml_Matrix_t* M) {
    cml_Matrix_t *result = malloc(sizeof(cml_Matrix_t));
    result->row_n = M->col_n;
    result->col_n = M->row_n;
    result->m = calloc(M->row_n*M->col_n, sizeof(double));
    int i, j;
    for(i=0;i<M->row_n;i++){
        for(j=0;j<M->col_n;j++){
            result->m[j*result->col_n+i] = M->m[cml_basicGetIndex(M, i, j)];
        }
    }
    return result;
}


cml_Matrix_t* cml_basicAdd(cml_Matrix_t *a, cml_Matrix_t *b) {
    int i, j;
    cml_Matrix_t *result = cml_basicZeros(a->row_n, a->col_n);
    for(i=0;i<a->row_n;i++){
        for(j=0;j<a->col_n;j++){
            int id = cml_basicGetIndex(result, i, j);
            result->m[id] = a->m[id] + b->m[id];
        }
    }
    return result;
}



cml_Matrix_t* cml_basicNumProd(double a, cml_Matrix_t* b) {
    cml_Matrix_t *result = cml_basicCopy(b);
    int i, j;
    for(i=0;i<b->row_n;i++){
        for(j=0;j<b->col_n;j++){
            int id = cml_basicGetIndex(b, i, j);
            result->m[id] = a * b->m[id];
        }
    }
    return result;
}


cml_Matrix_t* cml_basicMinus(cml_Matrix_t *a, cml_Matrix_t*b) {
    cml_Matrix_t *bb = cml_basicNumProd(-1., b);
    cml_Matrix_t *result = cml_basicAdd(a, bb);
    cml_basicFreeMatrix(&bb);
    return result;
}


cml_Matrix_t* cml_basicDot(cml_Matrix_t *a, cml_Matrix_t *b) {
    cml_Matrix_t * re = cml_basicZeros(a->row_n, b->col_n);
    int i, j, k;
    for(i=0;i<re->row_n;i++){
        for(j=0;j<re->col_n;j++){
            for(k=0;k<a->col_n;k++){
                re->m[cml_basicGetIndex(re,i,j)] += a->m[cml_basicGetIndex(a,i,k)] * b->m[cml_basicGetIndex(b,k,j)];
            }
        }
    }
    return re;
}


cml_Matrix_t* cml_basicSubMatrix(cml_Matrix_t* M, int row, int col) {
    cml_Matrix_t *re = cml_basicZeros(M->row_n-1, M->col_n-1);
    int i, j, k=0;
    for(i=0;i<M->row_n;i++){
        if(i == row)
            continue;
        for(j=0;j<M->col_n;j++){
            if(j == col)
                continue;
            re->m[k++] = M->m[cml_basicGetIndex(M, i, j)];
        }
    }
    return re;
}


double cml_basicDeterminant(cml_Matrix_t *a) {
    double result=0;
    if(a->row_n==1){
        result = a->m[0];
    }else if(a->row_n==2){
        result = a->m[0]*a->m[3] - a->m[1]*a->m[2];
    }else{
        for(int j=0;j<a->col_n;j++){
            if(a->m[cml_basicGetIndex(a,0,j)]==0){
                result += 0;
            }else{
                result += pow(-1, j+1+1)  * a->m[cml_basicGetIndex(a,0,j)] * cml_basicDeterminant(
                        cml_basicSubMatrix(a, 0, j));
            }
        }
    }
    return result;
}


cml_Matrix_t* cml_basicInverse(cml_Matrix_t* M) {
    if(M->row_n==1){
        return cml_basicNumProd(1./M->m[0], cml_basicOnesLike(M));
    }
    
    cml_Matrix_t *cofactor = cml_basicZerosLike(M);
    int i, j;
    for(i=0;i<cofactor->row_n;i++){
        for(j=0;j<cofactor->col_n;j++){
            cofactor->m[cml_basicGetIndex(cofactor,i,j)] = \
                pow(-1, i+j+2)*cml_basicDeterminant(cml_basicSubMatrix(M, i, j));
        }
    }
    return cml_basicNumProd(1./cml_basicDeterminant(M), cml_basicTranspose(cofactor));
}


cml_Matrix_t* cml_basicDiag(double* array, int l) {
    cml_Matrix_t *re = cml_basicZeros(l, l);
    for(int i=0;i<l;i++){
        re->m[cml_basicGetIndex(re, i, i)] = array[i];
    }
    return re;
}


cml_Matrix_t* cml_basicCross(cml_Matrix_t *a, cml_Matrix_t *b) {
    cml_Matrix_t *re = cml_basicZerosLike(a);
    re->m[0] = a->m[1]*b->m[2] - a->m[2]*b->m[1];
    re->m[1] = a->m[2]*b->m[0] - a->m[0]*b->m[2];
    re->m[2] = a->m[0]*b->m[1] - a->m[1]*b->m[0];
    return re;
}


cml_Matrix_t* cml_basicSlice(cml_Matrix_t* M, cml_Matrix_t* r, cml_Matrix_t* c) {
    cml_Matrix_t* re = cml_basicZeros(r->col_n, c->col_n);
    int i, j, p, q;
    for(i=0;i<r->col_n;i++){
        p = (int)r->m[i];
        for(j=0;j<c->col_n;j++){
            q = (int)c->m[j];
            re->m[cml_basicGetIndex(re, i, j)] = M->m[cml_basicGetIndex(M, p, q)];
        }
    }
    return re;
}


cml_Matrix_t* cml_basicRowSlice(cml_Matrix_t *M, cml_Matrix_t *RowSlice) {
    return cml_basicSlice(M, RowSlice, cml_basicRange(0, M->col_n, 1));
}


cml_Matrix_t* cml_basicColSlice(cml_Matrix_t *M, cml_Matrix_t *ColSlice) {
    return cml_basicSlice(M, cml_basicRange(0, M->row_n, 1), ColSlice);
}


cml_Matrix_t* cml_basicSum(cml_Matrix_t* M, int axis) {
    cml_Matrix_t *re = NULL;
    int i, j;
    if(axis==-1) {
        re = cml_basicZeros(1, 1);
        for(i=0;i<M->row_n*M->col_n;i++){
            re->m[0] += M->m[i];
        }
    }else if(axis==0) {
        re = cml_basicZeros(M->row_n, 1);
        for(i=0;i<M->row_n;i++){
            for(j=0;j<M->col_n;j++){
                re->m[i] += M->m[cml_basicGetIndex(M, i, j)];
            }
        }
    }else if(axis==1) {
        re = cml_basicZeros(1, M->col_n);
        for(j=0;j<M->col_n;j++){
            for(i=0;i<M->row_n;i++){
                re->m[j] += M->m[cml_basicGetIndex(M, i, j)];
            }
        }
    }else{
        re = NULL;
    }
    return re;
}


cml_Matrix_t* cml_basicMean(cml_Matrix_t* M, int axis) {
    cml_Matrix_t *re =  cml_basicSum(M, axis);
    if(axis==0){
        re = cml_basicNumProd(1/(double)M->col_n, re);
    }else if(axis==1){
        re = cml_basicNumProd(1/(double)M->row_n, re);
    }else if(axis==-1){
        re->m[0] /= (double)(M->row_n * M->col_n);
    }else{
        re = NULL;
    }
    return re;
}


cml_Matrix_t* cml_basicMatPow(cml_Matrix_t* M, int p) {
    // only square matrix can do this
    cml_Matrix_t *re = cml_basicIdentity(M->row_n);
    for(int i=0;i<p;i++){
        re = cml_basicDot(re, M);
    }
    return re;
}


cml_Matrix_t* cml_basicElemPow(cml_Matrix_t* M, double p) {
    cml_Matrix_t *re = cml_basicZerosLike(M);
    int i, j;
    for(i=0;i<M->row_n;i++){
        for(j=0;j<M->col_n;j++){
            re->m[cml_basicGetIndex(re, i, j)] = pow(M->m[cml_basicGetIndex(M, i, j)], p);
        }
    }
    return re;
}


cml_Matrix_t* cml_basicVecSort(cml_Matrix_t *vec) {
    cml_Matrix_t *re = cml_basicCopy(vec);
    int l = vec->row_n*(vec->row_n>=vec->col_n)
          + vec->col_n*(vec->col_n>vec->row_n);
    double temp;
    for(int i=0;i<l-1;i++){
        for(int j=i+1;j<l;j++){
            if(re->m[i]>re->m[j]){
                temp = re->m[i];
                re->m[i] = re->m[j];
                re->m[j] = temp;
            }
        }
    }
    return re;
}


cml_Matrix_t* cml_basicMatSort(cml_Matrix_t *M, int axis, int index) {
    cml_Matrix_t *vec = NULL;
    cml_Matrix_t *re = cml_basicCopy(M);
    int i, j, k;
    double tmp;
    if(axis==1){
        re = cml_basicTranspose(re);
    }
    vec = cml_basicRSlice(re, cml_basicRange(index, index+1, 1));
    for(j=0;j<re->col_n-1;j++){
        for(k=j;k<re->col_n;k++){
            if(vec->m[j]>vec->m[k]){
                for(i=0;i<re->row_n;i++){
                    tmp = re->m[cml_basicGetIndex(re,i,j)];
                    re->m[cml_basicGetIndex(re,i,j)] = re->m[cml_basicGetIndex(re,i,k)];
                    re->m[cml_basicGetIndex(re,i,k)] = tmp;
                }
            }
        }
    }
    if(axis==1){
        re = cml_basicTranspose(re);
    }
    cml_basicFreeMatrix(&vec);
    return re;
}


cml_Matrix_t* cml_basicReverse(cml_Matrix_t *vec) {
    cml_Matrix_t *re = cml_basicCopy(vec);
    int l = vec->row_n*(vec->row_n>=vec->col_n)
          + vec->col_n*(vec->col_n>vec->row_n);
    for(int i=0;i<l;i++){
        re->m[i] = vec->m[l-i-1];
    }
    return re;
}


double cml_basicVecMax(cml_Matrix_t *vec) {
    cml_Matrix_t *tmp = cml_basicVecSort(vec);
    double max;
    if(vec->row_n==1){
        max = tmp->m[vec->col_n-1];
    }else{
        max = tmp->m[vec->row_n-1];
    }
    cml_basicFreeMatrix(&tmp);
    return max;
}


double cml_basicVecMin(cml_Matrix_t *vec) {
    cml_Matrix_t *tmp = cml_basicVecSort(vec);
    double min = tmp->m[0];
    cml_basicFreeMatrix(&tmp);
    return min;
}


cml_Matrix_t* cml_basicMatMax(cml_Matrix_t *M, int axis) {
    cml_Matrix_t *re=NULL, *tmp=cml_basicCopy(M);

    if(axis==1){
        tmp = cml_basicTranspose(tmp);
    }
    re = cml_basicZeros(tmp->row_n, 1);
    
    for(int i=0;i<tmp->row_n;i++){
        re->m[i] = cml_basicVecMax(
            cml_basicRSlice(
                tmp, cml_basicRange(i,i+1,1)
            )
        );
        // printf("%d \n", re->m[i]);
    }
    if(axis==1){
        re = cml_basicTranspose(re);
    }else if(axis==-1){
        re->m[0] = cml_basicVecMax(re);
        re = cml_basicRSlice(re, cml_basicRange(0,1,1));
    }
    cml_basicFreeMatrix(&tmp);
    return re;
}


cml_Matrix_t* cml_basicMatMin(cml_Matrix_t *M, int axis) {
    cml_Matrix_t *re=NULL, *tmp=cml_basicCopy(M);

    if(axis==1){
        tmp = cml_basicTranspose(tmp);
    }
    re = cml_basicZeros(tmp->row_n, 1);
    
    for(int i=0;i<tmp->row_n;i++){
        re->m[i] = cml_basicVecMin(
            cml_basicRSlice(
                tmp, cml_basicRange(i,i+1,1)
            )
        );
    }
    if(axis==1){
        re = cml_basicTranspose(re);
    }
    if(axis==-1){
        re->m[0] = cml_basicVecMin(re);
        re = cml_basicRSlice(re, cml_basicRange(0,1,1));
    }
    cml_basicFreeMatrix(&tmp);
    return re;
}


int cml_basicApply(cml_Matrix_t* From, cml_Matrix_t* To, cml_Matrix_t* R_idx, cml_Matrix_t *C_idx) {
    // rp and cp should be range format and reformed
    int i, j, index;
    for(i=0;i<R_idx->col_n;i++){
        for(j=0;j<C_idx->col_n;j++){
            index = cml_basicGetIndex(To, R_idx->m[i], C_idx->m[j]);
            To->m[index] = From->m[cml_basicGetIndex(From, i, j)];
        }
    }
    return 1;
}


cml_Matrix_t* cml_basicShape(cml_Matrix_t *M) {
    cml_Matrix_t *re = cml_basicZeros(1, 2);
    re->m[0] = (double)M->row_n;
    re->m[1] = (double)M->col_n;
    return re;
}


cml_Matrix_t* cml_basicReshape(cml_Matrix_t *M, cml_Matrix_t *Shape){
    cml_Matrix_t *re = cml_basicCopy(M);
    re->row_n = (int)Shape->m[0];
    re->col_n = (int)Shape->m[1];
    return re;
}


cml_Matrix_t* cml_basicFlatten(cml_Matrix_t *M) {
    cml_Matrix_t *shape = cml_basicOnes(1, 2);
    shape->m[0] = 1;
    shape->m[1] = (double)M->row_n*M->col_n;
    cml_Matrix_t *re = cml_basicReshape(M, shape); 
    cml_basicFreeMatrix(&shape);
    return re;
}


cml_Matrix_t* cml_basicConcatenate(cml_Matrix_t *A, cml_Matrix_t *B, int axis) {
    int row=0, col=0;
    cml_Matrix_t *R_idx, *C_idx;
    cml_Matrix_t *new = cml_basicEmpty(1, A->row_n*A->col_n + B->row_n+B->col_n);

    if(axis==0) {
        row = A->row_n + B->row_n;
        col = A->col_n;
    }else if(axis==1) {
        col = A->col_n + B->col_n;
        row = A->row_n;
        A = cml_basicTranspose(A);
        B = cml_basicTranspose(B);
    }else {
        cml_basicFreeMatrix(&new);
        return NULL;
    }

    R_idx = cml_basicRange(0, A->row_n, 1);
    C_idx = cml_basicRange(0, B->col_n, 1);
    cml_basicApply(A, new, R_idx, C_idx);
    cml_basicFreeMatrix(&R_idx);
    cml_basicFreeMatrix(&C_idx);

    R_idx = cml_basicRange(A->row_n, A->row_n+B->row_n, 1);
    C_idx = cml_basicRange(A->col_n, A->col_n+B->col_n, 1);
    cml_basicApply(B, new, R_idx, C_idx);
    cml_basicFreeMatrix(&R_idx);
    cml_basicFreeMatrix(&C_idx);
    if(axis==1) {
        cml_Matrix_t *tmp = cml_basicTranspose(new);
        cml_basicFreeMatrix(&new);
        cml_basicFreeMatrix(&A);
        cml_basicFreeMatrix(&B);
        new = tmp;
    }
    return new;
}


#endif
