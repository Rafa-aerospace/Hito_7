from numpy import deg2rad, array, cos, sin, dot


def TBW_matrix(alpha: float, beta:float):

    """Creation of the rotation matrix Body Wind.

    Args:
        alpha (float): Angle of atack.          [Deg]
        beta (float): Sideslip.                 [Deg]

    Returns:
        Tbw(array): Rotation matrix.
    """

    alpha = deg2rad(alpha)
    beta = deg2rad(beta)

    Tbw = array([[cos(alpha)*cos(beta), -cos(alpha)*sin(beta), -sin(alpha)],
                    [sin(beta), cos(beta), 0],
                    [sin(alpha)*cos(beta), -sin(alpha)*sin(beta), cos(alpha)]])

    return Tbw

def TWB_matrix(alpha:float, beta:float):

    """Creation of the rotation matrix Wind Body.

    Args:
        alpha (float): Angle of atack.          [Deg]
        beta (float): Sideslip.                 [Deg]

    Returns:
        Tbw(array): Rotation matrix.
    """


    alpha = deg2rad(alpha)
    beta = deg2rad(beta)

    Twb = array([[cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta)],
                [-cos(alpha)*sin(beta), cos(beta), -sin(alpha)*cos(beta)],
                [-sin(alpha), 0 , cos(alpha)]])

    return Twb

# @njit(parallel=True)
def rotation_WB(alpha: float, beta: float, Data: array):
    """ Coordinates change from Wind frame to Body frame.

    Args:
        alpha (float): Angle of atack.          [Deg]
        beta (float): Sideslip.                 [Deg]
        Data (array): Data in the Wind frame.

    Returns:
        Rotated array: Data in the Body frame
    """

    Twb = TWB_matrix(alpha, beta)

    return dot(Twb,Data)

# @njit(parallel=True)
def rotation_BW(alpha: float, beta: float, Data: array):
    """ Coordinates change between Body and Wind frame.

    Args:
        alpha (float): Angle of atack           [Deg]
        beta (float): Sideslip                [Deg]
        Data (array): Data in the Body frame

    Returns:
        Rotated array: Data in the Wind frame
    """
    Tbw = TBW_matrix(alpha, beta)

    return dot(Tbw,Data)

def rotation_rev(phi: float, Data: array):
    """ Array rotation a certain angle of revolution around the X axis.

    Args:
        phi (float): Angle of revolution    [Deg]
        Data (array): Array to be rotated

    Returns:
        Rotated Array: Rotated array
    """
    phi = -deg2rad(phi) # QUÉ HACE UN MENOS AQUÍÍÍ AAAAAAAAAAAA

    Rot_x = array([[1 ,   0   ,    0    ],
                [0 , cos(phi), -sin(phi)],
                [0 , sin(phi),  cos(phi)]])

    return dot(Rot_x,Data)
