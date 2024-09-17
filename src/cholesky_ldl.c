//cholesky_ldl.c
/*
Autor: Emmanuel Alejandro Larralde Ortiz | emmanuel.larralde@cimat.mx
Descripcion:
    Factoriza una matriz cuadrada por el método Cholesky de la forma
        A = LDL', L' = transpose(L) y L[i][i] = 1

Compilar:
    gcc src/cholesky_ldl.c include/matrices/matrices.c -o output/cholesky_ldl.o
    chmod +x output/cholesky_ldl.o
Uso:
    ./output/cholesky_ldl.o <ruta_archivo_matriz> <ruta_archivo_txt> <verbose=0>
*/
#include <stdlib.h>
#include <stdio.h>
#include "../include/matrices/matrices.h"


int main(int argc, char **argv){
    //Cargando matriz
    double *mat;
    int m, n;
    mat = mat_from_txt(argv[1], &m, &n); //Lectura de la matriz
    if(mat == NULL){
        printf("No se pudo leer el archivo de matrices\n");
        return 0;
    }
    if(m != n){
        printf("No se puede factorizar matrices no cuadradas\n");
        return 0;
    }

    //Factorizacion ldl' = mat
    int n2 = n*n;
    double *l, *d;
    l = (double *) malloc(n2 * sizeof(double));
    d = (double *) malloc(n2 * sizeof(double));
    int check;
    check = cholesky_ldl(mat, l, d, n); //Funcion de la factorizacion
    if(check == -1){ //Method check
        printf("No es posible realizar una factorizacion LDL'\n");
        return 0;
    }

    //Comprobacion LDL' = A
    double *lt;
    lt = transpose(l, n, n); //L'
    double *aux_mat;
    aux_mat = matmul(l, d, n, n, n); //aux_mat = LD
    double *mat2;
    mat2 = matmul(aux_mat, lt, n, n, n); //mat2 = aux_matL' = LDL'
    double error;
    error = norm(mat, mat2, n2);
    printf("Error = norm(flatten(mat), flatten(LDL')): %lf\n", error);

    //Cargando vector b (Mx = b)
    double *b;
    int len;
    b = vec_from_txt(argv[2], &len); //Leer archivo del vector b
    if(len != n){
        printf("No se puede resolver sistema de ecuaciones con diferentes dimensiones\n");
        return 0;
    }

    //solve LDL'x = b
    //Solve LD(L'x) = b
    double *ltx;
    ltx = solve_l(aux_mat, b, n);
    //Solving L'x = ltx
    double *x;
    x = solve_u(lt, ltx, n);
    printf("x: ");
    print_matrix(x, 1, n);

    //Comprobacion LDL'x = b
    double *b2;
    b2 = matmul(mat2, x, n, n, 1);
    error = norm(b, b2, n);
    printf("Error = norm(b, b'): %lf\n", error);

    //Imprime resultados de la factorización si <verbose=1>
    if(argc > 3 && atoi(argv[3]) == 1){
        printf("L:\n");
        print_matrix(l, n, n);
        printf("D:\n");
        print_matrix(d, n, n);
    }

    free(mat);
    free(l);
    free(d);
    free(lt);
    free(aux_mat);
    free(mat2);
    free(b);
    free(ltx);
    free(x);
    free(b2);

    return 0;
}
