#Autor: Angel Octavio Parada Flores
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import math as mt 

def Fourier(B_0, A_m, B_m, m, k, y):
    VA = 0
    VB = 0
    for i in range(0, m):
        VA += A_m[i]*np.sin((i+1)*2*k*y)
        VB += B_m[i]*np.cos((i+1)*2*k*y)
    return B_0 + VA + VB 
    
x = sp.symbols("x")
m = int(input('Número de armónicos:', ))
I = float(input('Extremo inicial del periodo:',))
T = float(input('Periodo:', ))
F = sp.sympify(input('Función:', ))
p1 = input('¿Se va a reflejar la función? [S/n]:',)
p2 = 0
if p1 == 'n':
    p2 = input('¿Es una función en dos partes? [S/n]:',)
    
G = 0
if p1 == 'S':
    G = -F.subs(x, x - T/2)  
elif p1 == 'n' and p2 == 'n':
    G = F  
elif p1 == 'n' and p2 == 'S':
    G = sp.sympify(input('Función 2:', ))
    
k = mt.pi/T
B_0 = sp.integrate(F/T, (x, I, I+T/2)) + sp.integrate(G/T, (x, I+T/2, I+T))
A_m = []
B_m = []
for i in range (1, m+1):
    A_m.append(sp.integrate(F*(2/T)*sp.sin(i*2*k*x), (x, I, I+T/2)) + sp.integrate(G*(2/T)*sp.sin(i*2*k*x), (x, I+T/2, I+T)))
    B_m.append(sp.integrate(F*(2/T)*sp.cos(i*2*k*x), (x, I, I+T/2)) + sp.integrate(G*(2/T)*sp.cos(i*2*k*x), (x, I+T/2, I+T)))

y = np.arange(-T*3, T*3, 0.01)
plt.figure(figsize=(7,7))
plt.plot(y, Fourier(B_0, A_m, B_m, m, k, y), label = 'Núm. armónicos=' + repr(m))
plt.legend()
plt.grid()
plt.show()