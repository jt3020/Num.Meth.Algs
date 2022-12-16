import math 
import numpy as np 
import matplotlib.pyplot as plt

def Grid_Parabolic(xi,eta):
    r = 1 + xi 
    s = 1 + eta 
    x = (r**2 - s**2) / 2.0
    y = r * s
    return(x,y)

def Grid_Identity(xi,eta):
    x = xi 
    y = eta 
    return(x,y)

def Grid_Elliptic(xi,eta):
    a = 2
    r = 1 + xi
    s = np.pi * eta 
    x = a * np.cosh(r) * np.cos(s) 
    y = a * np.sinh(r) * np.sin(s)
    return(x,y)

def Grid_Horseshoe(xi,eta):
    rho = 2.0
    b0 =  1.0
    b1 =  2.0
    r = b0 + (b1 - b0) * eta 
    theta = np.pi * ( 1.0 - 2.0 * xi ) / 2.0
    x = rho * r * np.cos(theta) 
    y = r * np.sin(theta) 
    return(x,y)

Steps = 10

xi = np.zeros((Steps+1,Steps+1))
eta = np.zeros((Steps+1,Steps+1))
x = np.zeros((Steps+1,Steps+1))
y = np.zeros((Steps+1,Steps+1))

for i in range(0,Steps+1):
    for j in range(0,Steps+1):
        xi[i][j] = i*(1/Steps)
        eta[i][j] = j*(1/Steps)

# #Identity
# for i in range(0,Steps+1):
#     for j in range(0,Steps+1):      
#         x[i][j],y[i][j] = Grid_Identity(xi[i][j],eta[i][j])
#         plt.scatter(x[i][j],y[i][j],color='black',marker='x')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Identity')
# plt.show()

# #Parabolic
# for i in range(0,Steps+1):
#     for j in range(0,Steps+1):      
#         x[i][j],y[i][j] = Grid_Parabolic(xi[i][j],eta[i][j])
#         plt.scatter(x[i][j],y[i][j],color='black',marker='x')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Parabolic')
# plt.show()

#Elliptic
for i in range(0,Steps+1):
    for j in range(0,Steps+1):      
        x[i][j],y[i][j] = Grid_Elliptic(xi[i][j],eta[i][j])
        plt.scatter(x[i][j],y[i][j],color='black',marker='x')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Elliptic')
plt.show()

#Horseshoe
# for i in range(0,Steps+1):
#     for j in range(0,Steps+1):      
#         x[i][j],y[i][j] = Grid_Horseshoe(xi[i][j],eta[i][j])
#         plt.scatter(x[i][j],y[i][j],color='black',marker='x')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Horseshoe')
# plt.show()
