import numpy as np
import random

matrixSize = 400
A = np.random.rand(matrixSize, matrixSize)
B = np.dot(A, A.transpose())

print(matrixSize, matrixSize + 1)
for i in range(matrixSize):
    for j in range(matrixSize):
        print(str(round(B[j][i], 2)) + ' ', end='')
    print('1.0')
