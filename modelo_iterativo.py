"""Modelo termohidráulico de intercambiadores pillow-plate (PPHE).

Calcula, para varias geometrías de PPHE, el diámetro equivalente interno y
externo, las velocidades medias (obtenidas por el método iterativo de
iteracion_vel, ya que la ecuación (28) de O. Arsenyeva es implícita en w) y
los números adimensionales (Re, Pr, Nu) necesarios para obtener el factor de
fricción y el coeficiente de transferencia de calor por convección forzada
turbulenta.

Correlaciones tomadas de O. Arsenyeva et al.
"""

from dataclasses import dataclass
from math import sqrt
import formulas as form
import asignacion as asig
import geometría as geom
from iteracion_vel import resolver_velocidad
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
    U: float
    L_F: float


def imprimir_resultados(resultados):
    """Da formato a la impresión por consola de los resultados."""
    encabezado = (
        f"{'Caso':<8}{'Dh (in) [m]':>15}{'Dh (out) [m]':>15}{'u (in) [m/s]':>15}{'u (out) [m/s]':>15}{'Pr (in)':>10}{'Pr (out)':>10}{'Re (in)':>15}{'Re (out)':>15}{'Nu (in)':>12}{'Nu (out)':>12}{'f (in)':>12}{'f (out)':>12}{'h (in) [W/m2K]':>16}{'h (out) [W/m2K]':>16}{'U [W/m2K]':>14}{'L_F [m]':>12}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for res in resultados:
        print(
            f"{res.nombre:<8}{res.d_h_int:15.6e}{res.d_h_ext:15.6e}{res.u_m_int:15.4f}{res.u_m_ext:15.4f}{res.Pr_int:10.4f}{res.Pr_ext:10.4f}{res.Re_int:15.2e}{res.Re_ext:15.2e}{res.Nu_int:12.4f}{res.Nu_ext:12.4f}{res.f_int:12.4f}{res.f_ext:12.4f}{res.h_int:16.2f}{res.h_ext:16.2f}{res.U:14.2f}{res.L_F:12.4f}"
        )



# =============================================================================
# Cálculos numéricos
# =============================================================================

# Convenio de subíndices del solver (iteracion_vel): el fluido 1 es el que
# circula por el canal INTERNO, que es aquel cuya pérdida de carga se impone y
# cuya velocidad se despeja de la ec. (28). El fluido 2 (canal EXTERNO) es el
# CALIENTE: cede calor al interno, y su velocidad sale del balance de masa (30).

# Caso de estudio 1 del artículo (Tabla 4). Ojo con la numeración de esa tabla:
# su "fluido 1" es el CALIENTE, que aquí es el externo (ΔP° = 40 kPa), y su
# "fluido 2" es el FRÍO, que es el interno y el que impone la pérdida de carga
# de la ec. (28): ΔP° = 60 kPa. En la Tabla 5 el frío satura sus 60 kPa en las
# tres geometrías, mientras que el caliente se queda por debajo de sus 40 kPa.
              #    kg/m3       Pa*s       w/(m*K)   J/(kg*K)      °C         °C           Pa
props_fl_int = {"rho":980, "mu":0.0008, "k":0.618, "cp":4175, "t_in":10,  "t_out":50, "del_P":60e3}  #Fluido 2 del artículo: interno, frío
props_fl_ext = {"rho":980, "mu":0.0005, "k":0.654, "cp":4175, "t_in":70, "t_out":40}                #Fluido 1 del artículo: externo, caliente

resultados = []

for numero_iteracion, i in enumerate([geom.PPHE1, geom.PPHE2, geom.PPHE3], start=1):
    nombre = f"PPHE{numero_iteracion}"

    constantes_n = asig.asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i)

    """VELOCIDADES (método iterativo, ecuaciones (28) y (30) de O. Arsenyeva)"""

    print(f"\n--- {nombre} ---")
    u_m_int, u_m_ext, L_F, U_global = resolver_velocidad(
        rho1=props_fl_int.get("rho"), mu1=props_fl_int.get("mu"),
        cp1=props_fl_int.get("cp"), k1=props_fl_int.get("k"),
        n1_fric=constantes_n[0], n2_fric=constantes_n[1],
        n3_nu=constantes_n[2], n4_nu=constantes_n[3], n5_nu=constantes_n[4],
        rho2=props_fl_ext.get("rho"), mu2=props_fl_ext.get("mu"),
        cp2=props_fl_ext.get("cp"), k2=props_fl_ext.get("k"),
        t1in=props_fl_int.get("t_in"), t1out=props_fl_int.get("t_out"),
        t2in=props_fl_ext.get("t_in"), t2out=props_fl_ext.get("t_out"),
        del_P1=props_fl_int.get("del_P"),
        bi=i.b_i, b=i.b, wpp=i.w_pp, we=i.w_e, Fx=i.Fx, delta_w=i.delta_pp,
    )

    """CÁLCULO INTERNO"""

    d_h_int = form.deI(i.b_i)                                                                                                   # Diámetro equivalente
    Re_int = form.Re(rho=props_fl_int.get("rho"), u=u_m_int, dh=d_h_int, mu=props_fl_int.get("mu"))                              # Número de Reynolds
    Pr_int = form.Pr(cp=props_fl_int.get("cp"), mu=props_fl_int.get("mu"), k=props_fl_int.get("k"))                              # Número de Prandtl
    Nu_int = form.NuI(n3=constantes_n[2], n4=constantes_n[3], n5=constantes_n[4], Re=Re_int, Pr=Pr_int)                          # Número de Nusselt
    f_int = form.fI(n1=constantes_n[0], Re=Re_int, n2=constantes_n[1])                                                           # Factor de fricción
    h_int = Nu_int * props_fl_int.get("k") / d_h_int                                                                             # Coef. de transferencia de calor


    """CÁLCULO EXTERNO"""

    d_h_ext = form.deE(i.b_i, i.b)                                                                                              # Diámetro equivalente
    Re_ext = form.Re(rho=props_fl_ext.get("rho"), u=u_m_ext, dh=d_h_ext, mu=props_fl_ext.get("mu"))                              # Número de Reynolds
    Pr_ext = form.Pr(cp=props_fl_ext.get("cp"), mu=props_fl_ext.get("mu"), k=props_fl_ext.get("k"))                              # Número de Prandtl
    f_ext = form.fE(Re=Re_ext)                                                                                                  # Factor de fricción
    Nu_ext = form.NuE(Re=Re_ext, Pr=Pr_ext, f=f_ext)                                                                            # Número de Nusselt
    h_ext = Nu_ext * props_fl_ext.get("k") / d_h_ext                                                                             # Coef. de transferencia de calor

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
            U=U_global,
            L_F=L_F,
        )
    )

imprimir_resultados(resultados)
