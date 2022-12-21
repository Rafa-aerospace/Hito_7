
# g, Mach, gamma, alpha, beta
from Integrals import aero_coeffs
from Local_Inclination import Cp_3D, Cp_3D_Modified

def LIM_Application(g, Mach, gamma, alpha, beta):

    Normals_Tensor = g.N
    x = g.x
    z = g.z
    dx = g.dx
    dphi = g.dphi
    nx = g.nx
    nphi = g.nphi

    #Cp = Cp_3D(alpha, beta, nx, nphi, Normals_Tensor)
    Cp = Cp_3D_Modified(alpha, beta, nx, nphi, Normals_Tensor, Mach, gamma)

    CD_CY_CL = aero_coeffs(Cp, Normals_Tensor, x, z, nx, nphi, dphi, alpha, beta)
    

    return Cp , CD_CY_CL