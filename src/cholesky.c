//cholesky.c
/*
Autor: Emmanuel Alejandro Larralde Ortiz | emmanuel.larralde@cimat.mx
Descripcion:
    Factoriza una matriz cuadrada por el método Cholesky, es decir,
        A = LL', L' = transpose(L)
    El método de Cholesky está garantizado para matrices simétricas definidas positivas.

Compilar:
    gcc src/cholesky.c include/matrices/matrices.c -o output/cholesky.o
    chmod +x output/cholesky.o
Uso:
    ./output/cholesky.o <ruta_archivo_matriz> <ruta_archivo_txt> <verbose=0>
*/
#include <stdlib.h>
#include <stdio.h>
#include "../include/matrices/matrices.h"


int main(int argc, char **argv){
    double *mat;
    int m, n;
    mat = mat_from_txt(argv[1], &m, &n); //Leer matriz
    if(mat == NULL){
        printf("No se pudo leer el archivo de matrices\n");
        return 0;
    }
    if(m != n){
        printf("No se puede factorizar matrices no cuadradas\n");
        return 0;
    }

    //Factorizacion por metodo de cholesky
    int n2 = n*n;
    double *l;
    l = (double *) malloc(n2 * sizeof(double));
    int check;
    check = cholesky(mat, l, n); //factorizacion.
    if(check == -1){ //Method safeguard
        printf("No es posible realizar una factorizacion LL'\n");
        return 0;
    }

    //Comprobación LL' = mat
    double *lt;
    lt = transpose(l, n, n);
    double *mat2;
    mat2 = matmul(l, lt, n, n, n);
    double error;
    error = norm(mat, mat2, n2);
    printf("Error = norm(flatten(mat), flatten(LL')): %lf\n", error);

    //Cargando vector b (Mx = b)
    double *b;
    int len;
    b = vec_from_txt(argv[2], &len); //Lectura del vector b
    if(len != n){
        printf("No se puede resolver sistema de ecuaciones con diferentes dimensiones\n");
        return 0;
    }

    //Solving
    //LL'x = b
    //L(L'x) = b
    double *ltx;
    ltx = solve_l(l, b, n);

    //Solving L'x = ltx
    double *x;
    x = solve_u(lt, ltx, n);
    printf("x: ");
    print_matrix(x, 1, n);

    //Comprobando que Mx = b
    double *b2;
    b2 = matmul(mat2, x, n, n, 1);
    error = norm(b, b2, n);
    printf("Error = norm(b, b'): %lf\n", error);

    //Imprime resultados de la factorización si <verbose=1>
    if(argc > 3 && atoi(argv[3]) == 1){
        printf("L:\n");
        print_matrix(l, n, n);
    }

    free(mat);
    free(l);
    free(lt);
    free(mat2);
    free(ltx);
    free(x);
    free(b2);

    return 0;
}
