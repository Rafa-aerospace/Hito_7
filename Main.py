from Functionalities import Coefs_fM, Convergence, Basic


# "Horizontal_Parabola":
# "Semiesfera"
# "Cono"
# "F100"


nx = 200; nphi = 200
geometry_name = "F100"

Mach = 5; gamma = 1.4
alpha = 0; beta = 0

Basic(geometry_name, alpha, beta, Mach, nx, nphi,gamma)
# Coefs_fM(geometry_name, alpha, beta, nx, nphi)
# Convergence(geometry_name, Mach, alpha, beta)











