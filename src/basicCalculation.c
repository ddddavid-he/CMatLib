#ifndef BASIC_CALCULATION
#define BASIC_CALCULATION

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basicCalculation.h" 



void basicFreeMatrix(Matrix** Mp) {
    // CAUTION: Mp should be pointer of a Matrix Pointer.
    free((*Mp)->m);
    free(*Mp);
    *Mp = NULL;
}


int giveIndex(Matrix *M, int i, int j) {
    return M->col_n*i + j;
}


void basicShowMatrix(Matrix * M) {

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
            int index = giveIndex(M, i, j);
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


void basicFillMatrixFromArray(double* array, int row, int col, Matrix* M) {
    int i, j;
    M->row_n = row;
    M->col_n = col;
    double* m = (double*) calloc(row*col, sizeof(double));
    for(i=0;i<row;i++) {
        for(j=0;j<col;j++) {
            int index = giveIndex(M, i, j);
            m[index] = array[index];
        }
    }
    M->m = m;
}


Matrix* basicArray2Matrix(double* array, int row, int col) {
    Matrix *M = malloc(sizeof(Matrix));
    basicFillMatrixFromArray(array, row, col, M);
    return M;
}


Matrix* basicCopy(Matrix *M) {
    Matrix *result = malloc(sizeof(Matrix));
    result->m = calloc(M->row_n*M->col_n, sizeof(double));
    result->row_n = M->row_n;
    result->col_n = M->col_n;
    for(int i=0;i<M->row_n*M->col_n;i++){
        result->m[i] = M->m[i];
    }
    return result;
}


Matrix* basicRange(int from, int to, int skip) {
    // to is not included
    int length = (int)ceil(((double)(to-from))/((double)skip));
    Matrix * vec = basicZeros(1, length);
    if(skip>0){
        for(int i=from,j=0;i<to; i+=skip,j++){
            vec->m[j] = i;
        }
    }else{
        for(int i=from,j=0;i>to; i+=skip,j++){
            vec->m[j] = i;
        }
    }

    return vec;
} 


Matrix* basicZeros(int row, int col) {
    Matrix * M = malloc(sizeof(Matrix));
    M->row_n = row;
    M->col_n = col;
    M->m = calloc(row*col, sizeof(double));
    for(int i=0;i<row*col;i++) {
        M->m[i] = 0.;
    }
    return M;
}


Matrix* basicOnes(int row, int col) {
    Matrix * M = malloc(sizeof(Matrix));
    M->row_n = row;
    M->col_n = col;
    M->m = calloc(row*col, sizeof(double));
    for(int i=0;i<row*col;i++) {
        M->m[i] = 1.;
    }
    return M;
}


Matrix* basicZerosLike(Matrix *M) {
    return basicZeros(M->row_n, M->col_n);
}


Matrix* basicOnesLike(Matrix *M) {
    return basicOnes(M->row_n, M->col_n);
}


Matrix* basicIdentity(int d) {
    Matrix *re = basicZeros(d, d);
    for(int i=0;i<d;i++){
        re->m[giveIndex(re, i, i)] = 1.;
    }
    return re;
} 


Matrix* basicIdentityLike(Matrix* M) {
    // only square matrix can do this
    Matrix *re = basicZerosLike(M);
    for(int i=0;i<M->row_n;i++) {
        re->m[giveIndex(re, i, i)] = 1.;
    }
    return re;
}


Matrix* basicTranspose(Matrix* M) {
    Matrix *result = malloc(sizeof(Matrix));
    result->row_n = M->col_n;
    result->col_n = M->row_n;
    result->m = calloc(M->row_n*M->col_n, sizeof(double));
    int i, j;
    for(i=0;i<M->row_n;i++){
        for(j=0;j<M->col_n;j++){
            result->m[j*result->col_n+i] = M->m[giveIndex(M, i, j)];
        }
    }
    return result;
}


Matrix* basicAdd(Matrix *a, Matrix *b) {
    int i, j;
    Matrix *result = basicZeros(a->row_n, a->col_n);
    for(i=0;i<a->row_n;i++){
        for(j=0;j<a->col_n;j++){
            int id = giveIndex(result, i, j);
            result->m[id] = a->m[id] + b->m[id];
        }
    }
    return result;
}



Matrix* basicNumProd(double a, Matrix* b) {
    Matrix *result = basicCopy(b);
    int i, j;
    for(i=0;i<b->row_n;i++){
        for(j=0;j<b->col_n;j++){
            int id = giveIndex(b, i, j);
            result->m[id] = a * b->m[id];
        }
    }
    return result;
}


Matrix* basicMinus(Matrix *a, Matrix*b) {
    Matrix *bb = basicNumProd(-1., b);
    Matrix *result = basicAdd(a, bb);
    basicFreeMatrix(&bb);
    return result;
}


Matrix* basicDot(Matrix *a, Matrix *b) {
    Matrix * re = basicZeros(a->row_n, b->col_n);
    int i, j, k;
    for(i=0;i<re->row_n;i++){
        for(j=0;j<re->col_n;j++){
            for(k=0;k<a->col_n;k++){
                re->m[giveIndex(re,i,j)] += a->m[giveIndex(a,i,k)] * b->m[giveIndex(b,k,j)];
            }
        }
    }
    return re;
}


Matrix* __subMatrix__(Matrix* M, int r, int c) {
    Matrix *re = basicZeros(M->row_n-1, M->col_n-1);
    int i, j, k=0;
    for(i=0;i<M->row_n;i++){
        if(i==r)
            continue;
        for(j=0;j<M->col_n;j++){
            if(j==c)
                continue;
            re->m[k++] = M->m[giveIndex(M, i, j)];
        }
    }
    return re;
}


double basicDeterminant(Matrix *a) {
    double result=0;
    if(a->row_n==1){
        result = a->m[0];
    }else if(a->row_n==2){
        result = a->m[0]*a->m[3] - a->m[1]*a->m[2];
    }else{
        for(int j=0;j<a->col_n;j++){
            if(a->m[giveIndex(a,0,j)]==0){
                result += 0;
            }else{
                result += pow(-1, j+1+1)  * a->m[giveIndex(a,0,j)] * basicDeterminant(__subMatrix__(a, 0, j));
            }
        }
    }
    return result;
}


Matrix* basicInverse(Matrix* M) {
    if(M->row_n==1){
        return basicNumProd(1./M->m[0], basicOnesLike(M));
    }else{
        Matrix *cofactor = basicZerosLike(M);
    int i, j;
    for(i=0;i<cofactor->row_n;i++){
        for(j=0;j<cofactor->col_n;j++){
            cofactor->m[giveIndex(cofactor,i,j)] = \
                pow(-1, i+j+2)*basicDeterminant(__subMatrix__(M, i, j));
        }
    }
    return basicNumProd(1./basicDeterminant(M), basicTranspose(cofactor));
    }
}


Matrix* basicDiag(double* array, int l) {
    Matrix *re = basicZeros(l, l);
    for(int i=0;i<l;i++){
        re->m[giveIndex(re, i, i)] = array[i];
    }
    return re;
}


Matrix* basic_3d_cross(Matrix *a, Matrix *b) {
    Matrix *re = basicZerosLike(a);
    re->m[0] = a->m[1]*b->m[2] - a->m[2]*b->m[1];
    re->m[1] = a->m[2]*b->m[0] - a->m[0]*b->m[2];
    re->m[2] = a->m[0]*b->m[1] - a->m[1]*b->m[0];
    return re;
}


Matrix* basicSlice(Matrix* M, Matrix* r, Matrix* c) {
    Matrix* re = basicZeros(r->col_n, c->col_n);
    int i, j, p, q;
    for(i=0;i<r->col_n;i++){
        p = (int)r->m[i];
        for(j=0;j<c->col_n;j++){
            q = (int)c->m[j];
            re->m[giveIndex(re, i, j)] = M->m[giveIndex(M, p, q)];
        }
    }
    return re;
}


Matrix* basicRSlice(Matrix *M, Matrix *RowSlice) {
    return basicSlice(M, RowSlice, basicRange(0, M->col_n, 1));
}


Matrix* basicCSlice(Matrix *M, Matrix *ColSlice) {
    return basicSlice(M, basicRange(0, M->row_n, 1), ColSlice);
}


Matrix* basicSum(Matrix* M, int axis) {
    Matrix *re = NULL;
    int i, j;
    if(axis==-1) {
        re = basicZeros(1, 1);
        for(i=0;i<M->row_n*M->col_n;i++){
            re->m[0] += M->m[i];
        }
    }else if(axis==0) {
        re = basicZeros(M->row_n, 1);
        for(i=0;i<M->row_n;i++){
            for(j=0;j<M->col_n;j++){
                re->m[i] += M->m[giveIndex(M, i, j)];
            }
        }
    }else if(axis==1) {
        re = basicZeros(1, M->col_n);
        for(j=0;j<M->col_n;j++){
            for(i=0;i<M->row_n;i++){
                re->m[j] += M->m[giveIndex(M, i, j)];
            }
        }
    }else{
        re = NULL;
    }
    return re;
}


Matrix* basicMean(Matrix* M, int axis) {
    Matrix *re =  basicSum(M, axis);
    if(axis==0){
        re = basicNumProd(1/(double)M->col_n, re);
    }else if(axis==1){
        re = basicNumProd(1/(double)M->row_n, re);
    }else if(axis==-1){
        re->m[0] /= (double)(M->row_n * M->col_n);
    }else{
        re = NULL;
    }
    return re;
}


Matrix* basicMatPow(Matrix* M, int p) {
    // only square matrix can do this
    Matrix *re = basicIdentityLike(M);
    for(int i=0;i<p;i++){
        re = basicDot(re, M);
    }
    return re;
}


Matrix* basicElemPow(Matrix* M, double p) {
    Matrix *re = basicZerosLike(M);
    int i, j;
    for(i=0;i<M->row_n;i++){
        for(j=0;j<M->col_n;j++){
            re->m[giveIndex(re, i, j)] = pow(M->m[giveIndex(M, i, j)], p);
        }
    }
    return re;
}


Matrix* basicVecSort(Matrix *vec) {
    Matrix *re = basicCopy(vec);
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


Matrix* basicMatSort(Matrix *M, int axis, int index) {
    Matrix *vec = NULL;
    Matrix *re = basicCopy(M);
    int i, j, k;
    double tmp;
    if(axis==1){
        re = basicTranspose(re);
    }
    vec = basicRSlice(re, basicRange(index, index+1, 1));
    for(j=0;j<re->col_n-1;j++){
        for(k=j;k<re->col_n;k++){
            if(vec->m[j]>vec->m[k]){
                for(i=0;i<re->row_n;i++){
                    tmp = re->m[giveIndex(re,i,j)];
                    re->m[giveIndex(re,i,j)] = re->m[giveIndex(re,i,k)];
                    re->m[giveIndex(re,i,k)] = tmp;
                }
            }
        }
    }
    if(axis==1){
        re = basicTranspose(re);
    }
    basicFreeMatrix(&vec);
    return re;
}


Matrix* basicReverse(Matrix *vec) {
    Matrix *re = basicCopy(vec);
    int l = vec->row_n*(vec->row_n>=vec->col_n)
          + vec->col_n*(vec->col_n>vec->row_n);
    for(int i=0;i<l;i++){
        re->m[i] = vec->m[l-i-1];
    }
    return re;
}


double basicVecMax(Matrix *vec) {
    Matrix *tmp = basicVecSort(vec);
    double max;
    if(vec->row_n==1){
        max = tmp->m[vec->col_n-1];
    }else{
        max = tmp->m[vec->row_n-1];
    }
    basicFreeMatrix(&tmp);
    return max;
}


double basicVecMin(Matrix *vec) {
    Matrix *tmp = basicVecSort(vec);
    double min = tmp->m[0];
    basicFreeMatrix(&tmp);
    return min;
}


Matrix* basicMatMax(Matrix *M, int axis) {
    Matrix *re=NULL, *tmp=basicCopy(M);

    if(axis==1){
        tmp = basicTranspose(tmp);
    }
    re = basicZeros(tmp->row_n, 1);
    
    for(int i=0;i<tmp->row_n;i++){
        re->m[i] = basicVecMax(basicRSlice(tmp,basicRange(i,i+1,1)));
        // printf("%d \n", re->m[i]);
    }
    if(axis==1){
        re = basicTranspose(re);
    }else if(axis==-1){
        re->m[0] = basicVecMax(re);
        re = basicRSlice(re, basicRange(0,1,1));
    }
    basicFreeMatrix(&tmp);
    return re;
}


Matrix* basicMatMin(Matrix *M, int axis) {
    Matrix *re=NULL, *tmp=basicCopy(M);

    if(axis==1){
        tmp = basicTranspose(tmp);
    }
    re = basicZeros(tmp->row_n, 1);
    
    for(int i=0;i<tmp->row_n;i++){
        re->m[i] = basicVecMin(basicRSlice(tmp,basicRange(i,i+1,1)));
    }
    if(axis==1){
        re = basicTranspose(re);
    }else if(axis==-1){
        re->m[0] = basicVecMin(re);
        re = basicRSlice(re, basicRange(0,1,1));
    }
    basicFreeMatrix(&tmp);
    return re;
}


int basicApply(Matrix* From, Matrix* To, Matrix* rp, Matrix *cp) {
    // rp and cp should be range format and reformed
    int i, j, index;
    for(i=0;i<rp->col_n;i++){
        for(j=0;j<cp->col_n;j++){
            index = giveIndex(To, rp->m[i], cp->m[j]);
            To->m[index] = From->m[giveIndex(From, i, j)];
        }
    }
    return 1;
}


Matrix* basicShape(Matrix *M) {
    Matrix *re = basicZeros(1, 2);
    re->m[0] = (double)M->row_n;
    re->m[1] = (double)M->col_n;
    return re;
}


Matrix* basicReshape(Matrix *M, Matrix *Shape){
    Matrix *re = basicCopy(M);
    re->row_n = (int)Shape->m[0];
    re->col_n = (int)Shape->m[1];
    return re;
}


Matrix *basicFlatten(Matrix *M) {
    Matrix *shape = basicOnes(1, 2);
    shape->m[0] = 1;
    shape->m[1] = (double)M->row_n*M->col_n;
    Matrix *re = basicReshape(M, shape); 
    basicFreeMatrix(&shape);
    return re;
}



#endif
