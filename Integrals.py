from numpy import reshape, zeros, floor, sqrt, dot
from Rotations import rotation_BW

def aero_coeffs(Cp, N, x, z, n_x, n_phi, dphi, alpha, beta):

    N_panels = n_x*n_phi

    N_matrix = reshape(N, [N_panels, 3])
    Cp_v = reshape(Cp, [N_panels])

    d_Surface = zeros(N_panels)

    for n in range(N_panels):

        i = int(floor( n/ (n_phi) ))
        # print(i)

        hipotenusa = sqrt(( ( x[i+1] - x[i] )**2 + ( z[i+1] - z[i] )**2 ))

        d_Surface[n] = dphi * ( z[i+1] + z[i] )/2 * hipotenusa

    S_rev = sum(d_Surface)

    Coefs_XYZ = ( dot( ((-1) * Cp_v * d_Surface), N_matrix)  ) / S_rev

    CD_CY_CL =  rotation_BW(alpha, beta, Coefs_XYZ )

    return CD_CY_CL