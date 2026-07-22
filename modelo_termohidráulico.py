

def factor_de_fricción(constante1,reynods,constante2):
    """Esta función calcula el factor de fricción. Esta ecuación de ley de potencia se trata de una aproximación"""
    f = constante1 * reynods ** constante2
    return f

def nusselt(constante1,constante2,constante3, reynods, prandtl):
    """Esta función calcula el número de Nusselt. Esta ecuación de ley de potencia se trata de una aproximación"""
    nu = constante1 * reynods ** constante2 * prandtl ** constante3
    return nu

def reynolds(densidad, velocidad, diametro, viscosidad):
    """Cálculo del número de Reynolds"""
    reynolds = (densidad * velocidad * diametro) / viscosidad
    return reynolds

def prandtl(calor_especifico, viscosidad, conductividad):
    """Cálculo del número de Prandtl"""
    prandtl = (calor_especifico * viscosidad) / conductividad
    return prandtl

def asignacion_de_constantes(Re,sT,sL,dsp,h):
    if Re in range(1000,8000):
        a = 2*sL/sT
        b = dsp/sT
        c = h/sT
        if a in range(0.57,0.59) and b in range(0.1,0.14) and c in range(0.042,0.083):
            n1 = 8.74*b + (17*c + 0.73)
            n2 = -0.38
            n3 = 0.0775*b + (0.38*c + 0.005)
            n4 = 0.75
            n5 = 0.4
        elif a in range(0.99,1.01) and b in range(0.17,0.24) and c in range(0.071,0.143):
            n1 = -15.3*b + (1.4*c + 5.4)
            n2 = 1.725*b + (1.11*c - 0.66)
            n3 = 0.03*b + (0.76*c + 0.032)
            n4 = -1.12*c + 0.905
            n5 = 0.4
        elif a in range(1.7,1.72) and b in range(0.17,0.24) and c in range(0.071,0.143):
            n1 = 1.35*b + (2.8*c + 0.92)
            n2 = 0.3*b + (0.53*c - 0.29)
            n3 = -0.163*b + (0.711*c + 0.022)
            n4 = 0.29*b + (-c + 0.8)
            n5 = 0.4
        else:
            print("Las características geométricas no cumplen los requisitos")

    else:
        print("El número de Reynolds está fuera de rango para poder aplicar la correlación")

    return n1, n2, n3, n4, n5