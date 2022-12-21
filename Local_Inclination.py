from numpy import  array, dot, zeros, sin,  arccos, pi
from numpy.linalg import norm

from Rotations import rotation_BW
# from Normals import Rotated_Normals_Giant, Revoluting_Normals, Rotated_Normals_Loops

from OurDecorators import Function_Profiling
# from numba import njit
# @njit(parallel=True)
# def Theta_calculation(x, f):

#     theta = arctan( Analytical_Derivative( x, f ) )

#     theta[theta<0] = 0 # betta<0 implies Cp=0 based on Newton Theory

#     return theta
# def Cp_Dist_1D_by_Newton(x, f, M):
#     """Calculates the pressure coefficient over a 1D shape (defined by the evaluation of the function f over the parrtition x) at Mach M

#     Args:
#         x (array): linspace in the coordinate x
#         f (array): result of evaluating a function in x
#         M (float): Mach number

#     Returns:
#         _type_: _description_
#     """

#     Plot_geometry_1D(x, f(x))

#     Cp = 2*(sin(Theta_calculation(x, f)))**2

#     #Simple_plot_1D(x, Cp)

#     return Cp

# def Cp_Dist_1D_by_Newton_Modified(x, f, M, gamma=1.4, Mach=10):

#     Plot_geometry_1D(x, f(x), title=r'2D Geometry', axis_labels=[r'$x$ [-]', r'$z(x)$ [-]'])

#     Cp_max = 2/(gamma*M**2) * ( (((gamma+1)**2 * M**2)/(4*gamma*M**2-2*(gamma-1)))**(gamma/(gamma-1)) * (1-gamma+2*gamma*M**2)/(gamma+1) - 1 ) # From page 62 Andeson Hypersonic...

#     Cp_modified = Cp_max*( sin(Theta_calculation(x, f)) )**2

#     Simple_plot_1D(x, Cp_modified)

#     return Cp_modified
def Theta_calculation(Normal: array, V0: array):
    """Calculates the angle theta

    Args:
        Normal (array): Normal to the surface at that point
        V0 (array): Velocity vector represented in body frame

    Returns:
        Theta(float): angle [rad]
    """
    psi = arccos(dot(V0,Normal)/norm(V0))

    if (pi/2 - psi) < 0:
        theta = 0
    else:
        theta = pi/2 - psi

    return theta

@Function_Profiling
def Cp_3D(alpha: float, beta: float, nx:int, nphi: int, Normals_tensor):
    """This function calculates the Cp for each normal of the surface. Remember that the Normals_Tensor already contemplates both the angle of atack and the sideslip.

    Args:
        alpha (float): Angle of attack.       [Deg]
        beta (float): Sideslip.               [Deg]
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        Normals_tensor(array)

    Returns:

    Cp_matrix: Matrix (nx,nphi). Each row is a "x" partition and has has "nphi" Cps, this for every "x" partition (nx)
    """
    #Normals_tensor = Revoluting_Normals(nx, nphi, f)

    Cp_matrix = zeros([nx,nphi])

    V0 = array([-1,0,0])
    V0 = rotation_BW(alpha,beta,V0)

    for i in range(nx):
        for j in range(nphi):

            normal = (Normals_tensor[i,j])

            theta = Theta_calculation(normal, V0)

            Cp_matrix[i,j] = 2*(sin(theta))**2

    return Cp_matrix

@Function_Profiling
def Cp_3D_Modified(alpha: float, beta: float, nx:int, nphi: int, Normals_tensor, M: float, gamma: float):
    """This function calculates the modified Cp for each normal of the surface. Remember that the Normals_Tensor already contemplates both the angle of atack and the sideslip.

    Args:
        alpha (float): Angle of attack.       [Deg]
        beta (float): Sideslip.               [Deg]
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        Normals_tensor(array)
        M (float): : Mach number
        gamma (float): fluid's gamma

    Output:
        Cp_matrix: Matrix (nx,nphi). Each row is a "x" partition and has has "nphi" Cps, this for every "x" partition (nx)

    """
    Cp_matrix = zeros([nx,nphi])

    Cp_max = 2/(gamma*M**2) * ( (((gamma+1)**2 * M**2)/(4*gamma*M**2-2*(gamma-1)))**(gamma/(gamma-1)) * (1-gamma+2*gamma*M**2)/(gamma+1) - 1 )

    V0 = array([-1,0,0])
    V0 = rotation_BW(alpha,beta,V0)

    for i in range(nx):
        for j in range(nphi):

            normal = (Normals_tensor[i,j])

            theta = Theta_calculation(normal,V0)

            Cp_matrix[i,j] = Cp_max*(sin(theta))**2

    return  Cp_matrix



