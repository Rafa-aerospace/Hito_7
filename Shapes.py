# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 16:14:06 2022

@author: Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""

from numpy import logical_and, sqrt, hstack, linspace
# from OurDecorators import Plot_Simetric_Geometry_1D

def Shape_definition(name, nx):

    x = linspace(0, 1, nx+1)
    # dx = x[1]-x[0]

    if name == "Horizontal_Parabola":
        return Horizontal_Parabola(x), x
        # Horizontal_Parabola(x, title=name, x_label = "x", y_label = "y"), x
    elif name == "Semiesfera":
        return Quarter_of_a_Circle(x), x
    elif name == "Cono":
        return conito(x), x
    elif name == "F100":
        return Flechita(x), x
    else:
        print("Geometry has not been selected properly. Please select an available geometry.")


#### Funciones para definir la geometría ####
#@Plot_Simetric_Geometry_1D
def Horizontal_Parabola(x): # Libro Anderson, Cap 3, page 63

    return sqrt( (x+1)/0.769 )

# @Plot_Simetric_Geometry_1D
def Quarter_of_a_Circle(x):

    return sqrt(1 - (x-1)**2)

#@Plot_Simetric_Geometry_1D
def conito(x):

    return x

# @Plot_Simetric_Geometry_1D
def Flechita(x): # Se la ha inventado rafa

    f1 = x[x<=0.5]
    f2 = 2*x[logical_and(x>0.5, x<=0.75)]-0.5
    f3 = 1-4*(x[x>0.75]-0.75)

    return hstack((f1, f2, f3))
###########################