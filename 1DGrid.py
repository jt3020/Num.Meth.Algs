import math 
import numpy as np 
import matplotlib.pyplot as plt

def Grid_Compress(xi):
    L = 2
    x = (np.exp(L*xi)-1)/(np.exp(L)-1)
    return(x)

def Grid_Linear(xi,Start,End):
    x = (1-xi)*Start + xi*End
    return(x)

def Grid_Identity(xi):
    x = xi 
    return(x)

def Grid_Tangent(xi):
    x = np.tan(np.pi*xi/4.0)
    return(x)


Steps = 10

xi = np.zeros(Steps+1)
x = np.zeros(Steps+1)

for i in range(0,Steps+1):
    xi[i] = i*(1/Steps)


print('xi')
print(xi)

print('Linear x')
x = Grid_Linear(xi,0.0,1.0)
print(x)

print('Exponentail x')
x = Grid_Compress(xi)
print(x)

print('Identity x')
x = Grid_Identity(xi)
print(x)

print('Tangent x')
x = Grid_Tangent(xi)
print(x)

# plt.plot(x,xi)
# plt.show()
