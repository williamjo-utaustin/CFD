import numpy as np

L_x = 10
L_y = 10

N_x = 10
N_y = 10

h_x = L_x/N_x
h_y = L_y/N_y

x = np.zeros(N_x + 2)
y = np.zeros(N_y + 2)

for i in range(0, N_x + 2):
    x[i] = (i - 1) * h_x

for j in range(0, N_y + 2):
    y[j] = (j - 1) * h_y

print (x,y)