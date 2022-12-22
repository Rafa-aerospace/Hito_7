from Functionalities import Execute

Geometry_names = ["Horizontal_Parabola",
                  "Semiesfera",
                  "Cono",
                  "F100"]

Functionalities_names = ["Evolución de los Coeficientes con el Mach",
                         "Cáculo y representación 3-D",
                         "Convergencia con el número de elementos"]

#################################
nx = 200; nphi = 200
geometry_name = Geometry_names[1]

Mach = 5; gamma = 1.4
alpha = 0; beta = 0

application = Functionalities_names[1]
#################################

Execute(geometry_name, application, nx, nphi, alpha, beta, Mach, gamma)