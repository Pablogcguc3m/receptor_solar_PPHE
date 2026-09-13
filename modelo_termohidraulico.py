"""Modelo termohidráulico de intercambiadores pillow-plate (PPHE).

Calcula, para varias geometrías de PPHE, el diámetro hidráulico interno,
la velocidad media, y los números adimensionales (Re, Pr, Nu) necesarios
para obtener el factor de fricción y el coeficiente de transferencia de
calor por convección forzada turbulenta.

Correlaciones tomadas de O. Arsenyeva et al.
"""

from dataclasses import dataclass
from math import sqrt
import formulas as form
import asignacion as asig
from iteracion_vel import resolver_velocidad as w
# =============================================================================
# Estructuras de datos
# =============================================================================



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
# Cálculos numéricos
# =============================================================================

              #   kg/s     kg/m3     Pa*s    w/(m*K)  J/(kg*K)
props_fl_int = {"G":30, "rho":980, "mu":0.0008, "k":0.618, "cp":4175} #Propiedades del fluido interno
props_fl_ext = {"G":40, "rho":980, "mu":0.0005, "k":0.654, "cp":4175} #Propiedades del fluido externo

resultados = []

for numero_iteracion, i in enumerate([PPHE1, PPHE2, PPHE3], start=1):
    nombre = f"PPHE{numero_iteracion}"

    constantes_n = asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i)

    """CÁLCULO INTERNO"""

    G_in = props_fl_int.get("G")/i.n_pl                                                                                            #Gasto másico por placa

    d_h_int = diametro_hidraulico_interno(i.b_i)                                                                                # Diámetro hidráulico
    u_m_int = velocidad_media_int(G=G_in, rho=props_fl_int.get("rho"), b_i=i.b_i, w_pp=i.w_pp, w_e=i.w_e)              # Velocidad media
    Re_int = reynolds(rho=props_fl_int.get("rho"), u=u_m_int, dh=d_h_int, mu=props_fl_int.get("mu"))                                # Número de Reynolds
    Pr_int = prandtl(cp=props_fl_int.get("cp"), mu=props_fl_int.get("mu"), k=props_fl_int.get("k"))                                   # Número de Prandtl
    Nu_int = nusselt_int(n3=constantes_n[2], n4=constantes_n[3], n5=constantes_n[4], Re=Re_int, Pr=Pr_int)                      # Número de Nusselt
    f_int = factor_de_friccion_int(n1=constantes_n[0], Re=Re_int, n2=constantes_n[1])                                               # Factor de fricción
    h_int = Nu_int * props_fl_int.get("k") / d_h_int                                                                              # Coef. de transferencia de calor


    """CÁLCULO EXTERNO"""

    G_ex = props_fl_ext.get("G")/(i.n_pl-1)                                                                                            #Gasto másico por placa

    d_h_ext = diametro_hidraulico_externo(i.b_i,i.b)                                                                                # Diámetro hidráulico
    u_m_ext = velocidad_media_ext(G=G_ex, rho=props_fl_ext.get("rho"), b_i=i.b_i, b=i.b, w_pp=i.w_pp, w_e=i.w_e)              # Velocidad media
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
