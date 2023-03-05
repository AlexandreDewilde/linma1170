
def LU(A):
    n = len(A)
    for k in range(n):
        pivot = A[k][k]
        for i in range(k+1, n):
            A[i][k] /= pivot
            for j in range(k+1, n):
                A[i][j] -= A[i][k]*A[k][j]
            
def solve(A, b):
    LU(A)
    forward(A, b)
    backward(A, b)

def backward(A, b):
    n = len(b)
    for k in range(n-1, -1, -1):
        for i in range(k+1, n):
            b[k] -= A[k][i]*b[i]
        b[k] /= A[k][k]

def forward(A, b):
    n = len(A)
    for k in range(n):
        for i in range(k):
            b[k] -= A[k][i]*b[i]

p = [[777,915],[793,335]]
y = [386., 492.]
c = [[1,3,1], [2,-3,0], [1,1/3,0]]
# LU(p)
# print(p)
import numpy as np
print(np.linalg.solve(p, y))
# solve(p, y)
# print(y)