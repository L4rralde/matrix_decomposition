//lu_crout.c
/*
Autor: Emmanuel Alejandro Larralde Ortiz | emmanuel.larralde@cimat.mx
Descripción:
    Factoriza una matriz cuadrada en matrices LU y
    resuelve un sistema de ecuaciones de la forma:
        Mx = LUx = b

Compilar:
    gcc src/lu_crout.c include/matrices/matrices.c -o output/lu_crout.o
    chmod +x output/lu_crout.o

Uso:
    ./output/lu_crout.o <ruta_archivo_matriz> <ruta_archivo_txt> <verbose=0>
*/
#include <stdlib.h>
#include <stdio.h>
#include "../include/matrices/matrices.h"


int main(int argc, char **argv){
    double *mat;
    int m, n;

    mat = mat_from_txt(argv[1], &m, &n); //Leer matriz.

    if(mat == NULL){
        printf("No se pudo leer el archivo de matrices\n");
        return 0;
    }

    if(m != n){
        printf("No se puede factorizar matrices no cuadradas\n");
        return 0;
    }

    int n2 = n*n;
    double *l, *u;
    l = (double *) malloc(n2 * sizeof(double));
    u = (double *) malloc(n2 * sizeof(double));

    int check;
    check = lu_crout(mat, l, u, n);//Factorización lu
    if(check == - 1){
        printf("No es posible realizar una factorización LU\n");
        return 0;
    }

    //Comprobación LU = A
    double *mat2;
    mat2 = matmul(l, u, n, n, n); // mat2 = l*u
    double error;
    error = norm(mat, mat2, n*n);
    printf("Error = norm(flatten(mat), flatten(LU)); %lf\n", error);

    //Solving:
    //LUx = b
    //L(Ux) = b
    double *b;
    int len;
    b = vec_from_txt(argv[2], &len); //Leer vector b.
    if(len != n){
        printf("No se puede resolver sistema de ecuaciones con diferentes dimensiones\n");
        return 0;
    }

    double *ux;
    ux = solve_l(l, b, n);

    //Solving Ux = ux
    double *x;
    x = solve_u(u, ux, n);
    printf("x: ");
    print_matrix(x, 1, n);

    //Comprobación LUx = b
    double *b2;
    b2 = matmul(mat2, x, n, n, 1);
    error = norm(b, b2, n);
    printf("Error = norm(b, b'): %lf\n", error);

    //Imprime resultados de la factorización si <verbose=1>
    if(argc > 3 && atoi(argv[3]) == 1){
        printf("L:\n");
        print_matrix(l, n, n);
        printf("U:\n");
        print_matrix(u, n, n);
    }

    free(mat);
    free(mat2);
    free(l);
    free(u);
    free(b);
    free(x);
    free(b2);

    return 0;
}
