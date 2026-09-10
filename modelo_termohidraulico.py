"""Modelo termohidráulico de intercambiadores pillow-plate (PPHE).

Calcula, para varias geometrías de PPHE, el diámetro hidráulico interno,
la velocidad media, y los números adimensionales (Re, Pr, Nu) necesarios
para obtener el factor de fricción y el coeficiente de transferencia de
calor por convección forzada turbulenta.

Correlaciones tomadas de O. Arsenyeva et al.
"""

from dataclasses import dataclass
from math import sqrt


# =============================================================================
# Números adimensionales
# =============================================================================

def factor_de_friccion_int(n1, Re, n2):
    """Factor de fricción para el CANAL INTERNO (ley de potencia)."""
    return n1 * Re ** n2

def factor_de_friccion_ext(Re):
    """Factor de fricción para el CANAL EXTERNO (ley de potencia, se puede asumir cierta independientemente de la geometría con las constantes dadas)"""
    A = 2.187
    n = 0.356
    return A*Re**(-n)

def nusselt_int(n3, n4, n5, Re, Pr):
    """Número de Nusselt para el CANAL INTERNO (ley de potencia)."""
    return n3 * Re ** n4 * Pr ** n5

def nusselt_ext(Re, Pr, f):
    """Número de Nusselt para el CANAL EXTERNO (correlación experimental, se puede asumir cierta independientemente de la geometría)"""
    psi = 0.58
    return (psi*f/8*Re*Pr)/(1.07+12.7*(psi*f/8)**(1/2)*(Pr**(2/3)-1))

def reynolds(rho, u, dh, mu):
    """Número de Reynolds."""
    return (rho * u * dh) / mu


def prandtl(cp, mu, k):
    """Número de Prandtl."""
    return (cp * mu) / k


# =============================================================================
# Parámetros geométricos derivados
# =============================================================================

def velocidad_media_int(G, rho, b_i, w_pp, w_e):
    """Velocidad media por el conducto interno."""
    area_flujo = b_i / sqrt(2) * (w_pp - 2 * w_e)
    return 0.94 * (G / (rho * area_flujo))  # Se aplica un factor de corrección

def velocidad_media_ext(G, rho, b_i, b, w_pp, w_e):
    """Velocidad media por el conducto externo."""
    area_flujo = (b_i + b - b_i / sqrt(2)) * (w_pp - 2 * w_e)
    return 0.94 * (G / (rho * area_flujo))  # Se aplica un factor de corrección


def diametro_hidraulico_interno(b_i):
    """Diámetro hidráulico interno."""
    return 1.06 * 2 * (b_i / sqrt(2))  # Se aplica un factor de corrección

def diametro_hidraulico_externo(b_i, b):
    """Diámetro hidráulico externo."""
    return 2 * ((b_i+b)-b_i/sqrt(2))


# =============================================================================
# Constantes de correlación (n1..n5) según la geometría
# =============================================================================

def asignacion_de_constantes(sT, s2L, dsp, h):
    """Asigna valores a las constantes n1..n5 según los parámetros geométricos."""
    a = s2L / sT
    b = dsp / sT
    c = h / sT

    if 0.57 <= a <= 0.59 and 0.1 <= b <= 0.14 and 0.042 <= c <= 0.083:
        n1 = 8.74 * b + (17 * c + 0.73)
        n2 = -0.38
        n3 = 0.0775 * b + (0.38 * c + 0.005)
        n4 = 0.75
        n5 = 0.4
    elif 0.99 <= a <= 1.01 and 0.17 <= b <= 0.24 and 0.071 <= c <= 0.143:
        n1 = -15.3 * b + (1.4 * c + 5.4)
        n2 = 1.725 * b + (1.11 * c - 0.66)
        n3 = 0.03 * b + (0.76 * c - 0.032)
        n4 = -1.12 * c + 0.905
        n5 = 0.4
    elif 1.7 <= a <= 1.72 and 0.17 <= b <= 0.24 and 0.071 <= c <= 0.17:
        n1 = 1.35 * b + (2.8 * c + 0.92)
        n2 = 0.3 * b + (0.53 * c - 0.29)
        n3 = -0.163 * b + (0.711 * c + 0.022)
        n4 = 0.29 * b + (-c + 0.8)
        n5 = 0.4
    else:
        raise ValueError(
            f"Las características geométricas no cumplen los requisitos: a={a:.3f}, b={b:.3f}, c={c:.3f}"
        )

    return [n1, n2, n3, n4, n5]


# =============================================================================
# Estructuras de datos
# =============================================================================

@dataclass(frozen=True)
class PPHEGeometry:
    delta_pp: float  # Espesor de placa [m]
    b_i: float       # Altura/profundidad interna [m]
    b: float         # Distancia entre placas [m]
    s_2l: float      # Paso longitudinal [m]
    s_t: float       # Paso transversal [m]
    d_sp: float      # Diagonal de la soldadura [m]
    w_pp: float      # Anchura de la placa [m]
    l_pp: float      # Longitud de la placa [m]
    w_e: float       # Anchura de soldadura de borde [m]


@dataclass(frozen=True)
class PPHEResult:
    nombre: str
    d_h_int: float
    u_m_int: float
    Pr_int: float
    Re_int: float
    Nu_int: float
    f_int: float
    h_int: float
    d_h_ext: float
    u_m_ext: float
    Pr_ext: float
    Re_ext: float
    Nu_ext: float
    f_ext: float
    h_ext: float


def imprimir_resultados(resultados):
    """Da formato a la impresión por consola de los resultados."""
    encabezado = (
        f"{'Caso':<8}{'Dh (in) [m]':>15}{'Dh (out) [m]':>15}{'u (in) [m/s]':>15}{'u (out) [m/s]':>15}{'Pr (in)':>10}{'Pr (out)':>10}{'Re (in)':>15}{'Re (out)':>15}{'Nu (in)':>12}{'Nu (out)':>12}{'f (in)':>12}{'f (out)':>12}{'h (in) [W/m2K]':>16}{'h (out) [W/m2K]':>16}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for res in resultados:
        print(
            f"{res.nombre:<8}{res.d_h_int:15.6e}{res.d_h_ext:15.6e}{res.u_m_int:15.4f}{res.u_m_ext:15.4f}{res.Pr_int:10.4f}{res.Pr_ext:10.4f}{res.Re_int:15.2e}{res.Re_ext:15.2e}{res.Nu_int:12.4f}{res.Nu_ext:12.4f}{res.f_int:12.4f}{res.f_ext:12.4f}{res.h_int:16.2f}{res.h_ext:16.2f}"
        )


# =============================================================================
# Geometrías a estudiar (Tabla 1, artículo de O. Arsenyeva et al.)
# =============================================================================

PPHE1 = PPHEGeometry(
    delta_pp=0.8e-3,
    b_i=3.4e-3,
    b = 2*3.4e-3,
    s_2l=42e-3,
    s_t=72e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3,
)
PPHE2 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=3.0e-3,
    b = 2*3.0e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3,
)
PPHE3 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=7.0e-3,
    b = 2*7.0e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3,
)

# =============================================================================
# Cálculos numéricos
# =============================================================================

              #   kg/s     kg/m3     Pa*s    w/(m*K)  J/(kg*K)
props_fl_int = {"G":30/70, "rho":980, "mu":0.0008, "k":0.618, "cp":4175} #Propiedades del fluido interno
props_fl_ext = {"G":40/69, "rho":980, "mu":0.0005, "k":0.654, "cp":4175} #Propiedades del fluido externo

resultados = []

for numero_iteracion, i in enumerate([PPHE1, PPHE2, PPHE3], start=1):
    nombre = f"PPHE{numero_iteracion}"

    constantes_n = asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i)

    """CÁLCULO INTERNO"""

    d_h_int = diametro_hidraulico_interno(i.b_i)                                                                                # Diámetro hidráulico
    u_m_int = velocidad_media_int(G=props_fl_int.get("G"), rho=props_fl_int.get("rho"), b_i=i.b_i, w_pp=i.w_pp, w_e=i.w_e)              # Velocidad media
    Re_int = reynolds(rho=props_fl_int.get("rho"), u=u_m_int, dh=d_h_int, mu=props_fl_int.get("mu"))                                # Número de Reynolds
    Pr_int = prandtl(cp=props_fl_int.get("cp"), mu=props_fl_int.get("mu"), k=props_fl_int.get("k"))                                   # Número de Prandtl
    Nu_int = nusselt_int(n3=constantes_n[2], n4=constantes_n[3], n5=constantes_n[4], Re=Re_int, Pr=Pr_int)                      # Número de Nusselt
    f_int = factor_de_friccion_int(n1=constantes_n[0], Re=Re_int, n2=constantes_n[1])                                               # Factor de fricción
    h_int = Nu_int * props_fl_int.get("k") / d_h_int                                                                              # Coef. de transferencia de calor


    """CÁLCULO EXTERNO"""

    d_h_ext = diametro_hidraulico_externo(i.b_i,i.b)                                                                                # Diámetro hidráulico
    u_m_ext = velocidad_media_ext(G=props_fl_ext.get("G"), rho=props_fl_ext.get("rho"), b_i=i.b_i, b=i.b, w_pp=i.w_pp, w_e=i.w_e)              # Velocidad media
    Re_ext = reynolds(rho=props_fl_ext.get("rho"), u=u_m_ext, dh=d_h_ext, mu=props_fl_ext.get("mu"))                                # Número de Reynolds
    Pr_ext = prandtl(cp=props_fl_ext.get("cp"), mu=props_fl_ext.get("mu"), k=props_fl_ext.get("k"))                                   # Número de Prandtl
    f_ext = factor_de_friccion_ext(Re=Re_ext)                                               # Factor de fricción
    Nu_ext = nusselt_ext(Re=Re_ext, Pr=Pr_ext, f=f_ext)                      # Número de Nusselt
    h_ext = Nu_ext * props_fl_ext.get("k") / d_h_ext                                                                              # Coef. de transferencia de calor

    resultados.append(
        PPHEResult(
            nombre=nombre,
            d_h_int=d_h_int,
            u_m_int=u_m_int,
            Pr_int=Pr_int,
            Re_int=Re_int,
            Nu_int=Nu_int,
            f_int=f_int,
            h_int=h_int,
            d_h_ext=d_h_ext,
            u_m_ext=u_m_ext,
            Pr_ext=Pr_ext,
            Re_ext=Re_ext,
            Nu_ext=Nu_ext,
            f_ext=f_ext,
            h_ext=h_ext,
        )
    )

imprimir_resultados(resultados)
