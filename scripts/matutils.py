#!/usr/bin/env python3
"""
Autor: Emmanuel Alejandro Larralde Ortiz | emmanuel.larralde@cimat.mx
Descripcion:
    Modulo con varias funciones para leer archivos de matrices y vectores,
    y funciones miscelaneas.
"""
import numpy as np

def read_matrix(fname: str) -> np.array:
    """
    Lee un archivo de una matriz.
    """
    with open(fname, 'r') as f:
        lines = f.read().splitlines()
    n,m = lines[0].split(" ")
    matrix = [[0.0 for _ in range(int(m))] for _ in range(int(n))]
    for i in range(int(n)):
        row = lines[i+1].split(' ')
        for j in range(int(m)):
            matrix[i][j] = float(row[j])
    return np.array(matrix)

def read_vector(fname: str) -> np.array:
    """
    Lee un archivo de un vector.
    """
    with open(fname, 'r') as f:
        lines = f.read().splitlines()
    n = int(lines.pop(0))
    return np.array([float(x) for x in lines[-1].split(" ")])

def write_matrix(matrix: np.array, fname: str) -> None:
    """
    Escribe un arreglo de numpy en un archivo de texto
    """
    shape_str = ' '.join([str(n) for n in matrix.shape])
    try:
        mat_str = [' '.join([str(x) for x in row]) for row in matrix]
    except:
        mat_str = [' '.join([str(x) for x in matrix])]
    with open(fname, 'w') as f:
        f.write(f"{shape_str}\n")
        for row_str in mat_str:
            f.write(f"{row_str}\n")


def positive_semidefinite(n: int) -> np.array:
    """
    Genera una matriz semi-definida positiva de tamaño n x n.
    Por definicion, el prioducto AA' es semi-deifinida positiva.
    A' = transpose(A)
    """
    A = np.random.rand(n, n)
    return np.dot(A, A.transpose())
