from dataclasses import dataclass
from math import sqrt

### Definición de funciones para calcular números adimensionales, parámetros para cálculo de estos y constantes de correlación

def factor_de_fricción(n1,Re,n2):
    """Esta función calcula el factor de fricción. Esta ecuación de ley de potencia se trata de una aproximación"""
    return n1 * Re ** n2

def nusselt(n3,n4,n5, Re, Pr):
    """Esta función calcula el número de Nusselt. Esta ecuación de ley de potencia se trata de una aproximación"""
    return n3 * Re ** n4 * Pr ** n5


def reynolds(rho, u, dh, mu):
    """Cálculo del número de Reynolds"""
    return (rho * u * dh) / mu

def prandtl(cp, mu, k):
    """Cálculo del número de Prandtl"""
    return (cp * mu) / k

def velocidad_media(G, rho, b_i, w_pp, w_e):
    """Cálculo de la velocidad media por el conducto interno"""
    area_flujo = b_i / sqrt(2) * (w_pp - 2 * w_e)
    return 0.94 * (G / (rho * area_flujo))  # Se aplica un factor de corrección

def diametro_hidraulico_interno(b_i):
    """Cálculo del diámetro hidráulico interno"""
    return 1.06 * 2 * (b_i/sqrt(2)) # Se aplica un factor de corrección

def asignacion_de_constantes(sT,s2L,dsp,h):
    """Esta función se encarga de asginar valores a las constantes ni"""
    a = s2L/sT
    b = dsp/sT
    c = h/sT
    if 0.57 <= a <= 0.59 and 0.1 <= b <= 0.14 and 0.042 <= c <= 0.083:
        n1 = 8.74*b + (17*c + 0.73)
        n2 = -0.38
        n3 = 0.0775*b + (0.38*c + 0.005)
        n4 = 0.75
        n5 = 0.4
    elif 0.99 <= a <= 1.01 and 0.17 <= b <= 0.24 and 0.071 <= c <= 0.143:
        n1 = -15.3*b + (1.4*c + 5.4)
        n2 = 1.725*b + (1.11*c - 0.66)
        n3 = 0.03*b + (0.76*c + 0.032)
        n4 = -1.12*c + 0.905
        n5 = 0.4
    elif 1.7 <= a <= 1.72 and 0.17 <= b <= 0.24 and 0.071 <= c <= 0.17:
        n1 = 1.35*b + (2.8*c + 0.92)
        n2 = 0.3*b + (0.53*c - 0.29)
        n3 = -0.163*b + (0.711*c + 0.022)
        n4 = 0.29*b + (-c + 0.8)
        n5 = 0.4
    else:
        raise ValueError(
            f"Las características geométricas no cumplen los requisitos: a={a:.3f}, b={b:.3f}, c={c:.3f}"
        )

    return [n1, n2, n3, n4, n5]

#### Parámetros geométricos para calcular. Sacados de la Tabla 1 del artículo de o. Arsenyeva et al.

@dataclass(frozen=True)
class PPHEGeometry:
    delta_pp: float   # Espesor de placa [m]
    b_i: float        # Altura/profundidad interna [m]
    s_2l: float       # Paso longitudinal [m]
    s_t: float        # Paso transversal [m]
    d_sp: float        # Diagonal de la soldadura [m]
    w_pp: float       # Anchura de la placa [m]
    l_pp: float       # Longitud de la placa [m]
    w_e: float        # Anchura de soldadura de borde [m]

@dataclass(frozen=True)
class PPHEResult:
    nombre: str
    d_h: float
    u_m: float
    Pr: float
    Re: float
    Nu: float
    f: float
    h: float


def imprimir_resultados(resultados):
    encabezado = (
        f"{'Caso':<8}{'Dh [m]':>12}{'u [m/s]':>12}{'Pr':>10}{'Re':>15}{'Nu':>12}{'f':>12}{'h [W/m2K]':>16}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for res in resultados:
        print(
            f"{res.nombre:<8}{res.d_h:12.6e}{res.u_m:12.4f}{res.Pr:10.4f}{res.Re:15.2e}{res.Nu:12.4f}{res.f:12.4f}{res.h:16.2f}"
        )

PPHE1 = PPHEGeometry(
    delta_pp=0.8e-3,
    b_i=3.4e-3,
    s_2l=42e-3,
    s_t=72e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3)
PPHE2 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=3.0e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3)
PPHE3 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=7.0e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3)

### CÁLCULOS NUMÉRICOS
G = 30 #kg/s
rho = 980 #kg/m3
mu = 0.0008 #Pa*s
k = 0.618 #W/(m*K)
cp = 4175 #J/(kg*K)
resultados = []

for numero_iteracion, i in enumerate([PPHE1, PPHE2, PPHE3], start=1):
    nombre = f"PPHE{numero_iteracion}"
    constantes_n = asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i)
    d_h = diametro_hidraulico_interno(i.b_i)  # Cálculo del diámetro hidráulico
    u_m = velocidad_media(G=G, rho=rho, b_i=i.b_i, w_pp=i.w_pp, w_e=i.w_e)  # Cálculo de la velocidad media
    Re = reynolds(rho=rho, u=u_m, dh=d_h, mu=mu)  # Cálculo del número de Reynolds
    Pr = prandtl(cp=cp, mu=mu, k=k)  # Cálculo del número de Prandtl
    Nu = nusselt(n3=constantes_n[2], n4=constantes_n[3], n5=constantes_n[4], Re=Re, Pr=Pr)  # Cálculo del número de Nusselt
    f = factor_de_fricción(n1=constantes_n[3], Re=Re, n2=constantes_n[4])  # Cálculo del factor de fricción
    h = Nu * k / d_h  # Cálculo del coeficiente de transferencia de calor
    resultados.append(
        PPHEResult(
            nombre=nombre,
            d_h=d_h,
            u_m=u_m,
            Pr=Pr,
            Re=Re,
            Nu=Nu,
            f=f,
            h=h,
        )
    )

imprimir_resultados(resultados)