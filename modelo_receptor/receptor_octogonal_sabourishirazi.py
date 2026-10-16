"""Receptor octogonal / en estrella con el Nusselt interno de Sabourishirazi.

Hermano de receptor_octogonal.py: mismo receptor, mismo mapa, mismo circuito
térmico, misma numérica y misma pérdida de carga (factor de fricción de Piper,
ec. 18 de Arsenyeva con n1, n2 de asignacion). Lo ÚNICO que cambia es el número
de Nusselt del canal interno, que pasa de la ec. (20) con las constantes de
Piper a la ec. (8) de

    A. Sabourishirazi, J.-L. Liow, M. Ghodrat, "Estimations of heat transfer
    coefficients in the channels of pillow-plate heat exchangers", Int. J. Heat
    Mass Transf. 251 (2025) 127421.

    Nu = a1 * (d*_ic,max)^a2 * (d*_ws)^a4 * (theta**)^a5 * (s_R)^a6 * Re^b * Pr^0.4

con los parámetros adimensionales de su sección 4:

    s_D       = sqrt((s_T/2)^2 + s_L^2)              paso diagonal
    s_D,eff   = s_D - d_ws                           su tramo efectivo
    theta     = (180/pi)*atan(s_T/(2 s_L))           ángulo del paso diagonal [grados]
    d*_ic,max = delta_ic,max / s_D,eff               (delta_ic,max = b_i)
    d*_ws     = d_ws / s_D,eff                       (d_ws = d_sp)
    theta**   = theta / 90
    s_R       = (2 s_L - d_ws) / (s_T - d_ws)

TIPO DE CANAL. El artículo clasifica por los pasos: 2s_L > s_T longitudinal,
2s_L < s_T transversal, 2s_L = s_T equidistante (su Tabla 4 tiene una fila de
coeficientes para cada uno). OJO, no coincide con los nombres de asignacion.py:
la "celda longitudinal" de allí (2s_L/s_T = 0.58, la del PPHE1) es un canal
TRANSVERSAL para Sabourishirazi, y la "celda transversal" (1.71, PPHE2/3) es
LONGITUDINAL. Se usan las tres filas de la Tabla 4 y se elige la del tipo de la
geometría; en el equidistante theta** = 0.5 y s_R = 1 son constantes y el
artículo los quita (a5 = a6 = 0).

REYNOLDS Y DIÁMETRO. El artículo define Re y Nu con el diámetro hidráulico
medio 4V/A_wet de su simulación de hidroconformado, que no tiene expresión
cerrada aquí. Se mantienen la velocidad de la ec. (3), la sección de la (16) y
el d_e de la (14), exactamente como en receptor_octogonal.py, para que la
comparación entre los dos métodos sólo vea el cambio de correlación.

VALIDEZ (Tablas 4 y 6). Pr entre 1 y 13; Re de 1200 a 20500 en el longitudinal
y de 1000 a 10500 en los otros dos; y los cuatro adimensionales dentro de los
rangos de la Tabla 6. Fuera de ellos se calcula igual, pero se avisa. Con
re_limitado = True, el Nu se evalúa con Re recortado al máximo de la Tabla 4 de
su tipo de canal, y la fricción de Piper con Re recortado a 8000 (ver
receptor_octogonal.comparacion_rango).
"""

from dataclasses import dataclass, field, replace
from math import atan, degrees, pi, sqrt
from typing import NamedTuple

import numpy as np

import receptor_octogonal as ro
from receptor_octogonal import (ALTURA, BARRIDO_LAMBDA, C, LAMBDA, LAMBDAS, N_PLACAS,
                                Q_MEDIO, Q_PICO, MapaMitad, ancho_placa)


# =============================================================================
# Correlación de Sabourishirazi (ec. 8, Tablas 4 y 6)
# =============================================================================

PR_EXPONENTE = 0.4
PR_RANGO = (1.0, 13.0)

# tipo: (a1, a2, a4, a5, a6, b), rango de Re, rangos de la Tabla 6
# (d*_ic,max, d*_ws, theta**, s_R); None donde el término no interviene
TABLA_4 = {
    "longitudinal": dict(
        a=(0.0465, 0.3765, -0.0898, 0.3989, 0.5984, 0.8549), Re=(1200.0, 20500.0),
        rangos=((0.0869, 0.2702), (0.2182, 0.625), (0.2473, 0.4229), (1.3240, 3.9173))),
    "transversal": dict(
        a=(0.3947, 0.3867, 0.1339, 1.9544, 0.9653, 0.7703), Re=(1000.0, 10500.0),
        rangos=((0.0869, 0.2702), (0.2182, 0.625), (0.5771, 0.7527), (0.2553, 0.7553))),
    "equidistante": dict(
        a=(0.0750, 0.3078, 0.1737, 0.0, 0.0, 0.8052), Re=(1000.0, 10500.0),
        rangos=((0.0749, 0.3447), (0.1071, 1.0), None, None)),
}
TOL_EQUIDISTANTE = 0.01   # |2s_L/s_T - 1| por debajo de esto, equidistante


class Sabourishirazi(NamedTuple):
    """Parámetros adimensionales de la geometría y la fila de la Tabla 4 que toca."""
    tipo: str
    s_D_eff: float        # [m]
    theta: float          # [grados]
    d_ic: float           # d*_ic,max
    d_ws: float           # d*_ws
    theta_2: float        # theta**
    s_R: float
    a: tuple              # (a1, a2, a4, a5, a6, b)
    Re: tuple             # Rango de Re de la Tabla 4
    fuera: tuple          # Nombres de los adimensionales fuera de la Tabla 6

    @property
    def factor(self):
        """Todo lo que no es Re^b*Pr^0.4 en la ec. (8)."""
        a1, a2, a4, a5, a6, _ = self.a
        return (a1 * self.d_ic ** a2 * self.d_ws ** a4
                * self.theta_2 ** a5 * self.s_R ** a6)


def parametros_sabourishirazi(s_t, s_2l, d_ws, delta_ic):
    """Adimensionales de la sección 4 del artículo y tipo de canal. Todo en metros."""
    s_l = s_2l / 2
    s_D_eff = sqrt((s_t / 2) ** 2 + s_l ** 2) - d_ws
    if s_D_eff <= 0:
        raise ValueError("El punto de soldadura no cabe en el paso diagonal (s_D <= d_ws).")
    theta = degrees(atan(s_t / s_2l))
    relacion = s_2l / s_t
    if abs(relacion - 1.0) <= TOL_EQUIDISTANTE:
        tipo = "equidistante"
    else:
        tipo = "longitudinal" if relacion > 1.0 else "transversal"
    fila = TABLA_4[tipo]
    valores = (delta_ic / s_D_eff, d_ws / s_D_eff, theta / 90.0,
               (s_2l - d_ws) / (s_t - d_ws))
    nombres = ("d*_ic,max", "d*_ws", "theta**", "s_R")
    fuera = tuple(nom for nom, v, r in zip(nombres, valores, fila["rangos"])
                  if r is not None and not r[0] <= v <= r[1])
    return Sabourishirazi(tipo, s_D_eff, theta, *valores, a=fila["a"], Re=fila["Re"],
                          fuera=fuera)


def NuI_sabourishirazi(s, Re, Pr):
    """Nusselt del canal interno, ec. (8) de Sabourishirazi. Admite arrays."""
    return s.factor * Re ** s.a[5] * Pr ** PR_EXPONENTE


# =============================================================================
# El receptor: el de receptor_octogonal.py con otro h_interno
# =============================================================================


@dataclass(frozen=True)
class ReceptorSabourishirazi(ro.ReceptorOctogonal):
    """ReceptorOctogonal con el Nusselt interno de Sabourishirazi (ec. 8)."""

    sab: Sabourishirazi = field(init=False)

    def __post_init__(self):
        super().__post_init__()
        i = self.panel
        object.__setattr__(self, "sab", parametros_sabourishirazi(
            s_t=i.s_t, s_2l=i.s_2l, d_ws=i.d_sp, delta_ic=i.b_i))

    nombre_nu = "Sabourishirazi et al. (2025), ec. (8)"

    @property
    def rango_re_nu(self):
        """Rango de Re de la Tabla 4 para el tipo de canal de la geometría."""
        return self.sab.Re

    def nusselt(self, Re, Pr):
        """Nusselt del canal interno, ec. (8). Lo único que cambia respecto de Piper."""
        return NuI_sabourishirazi(self.sab, Re, Pr)


# =============================================================================
# Salida por consola
# =============================================================================


def validez(res):
    """Re y Pr de todas las porciones frente a los rangos de las Tablas 4 y 6."""
    r, s = res.receptor, res.receptor.sab
    Re, Pr = [], []
    for p in res.placas:
        # Igual que en la pérdida de carga: cada porción a su temperatura media
        entradas = np.empty_like(p.T_fluido)
        orden = range(r.N) if p.sube else range(r.N - 1, -1, -1)
        T = np.full(r.M, p.T_ent)
        for j in orden:
            entradas[j], T = T, p.T_fluido[j]
        pr = r.fluido.propiedades(0.5 * (entradas + p.T_fluido))
        Re.append(r.h_interno(0.5 * (entradas + p.T_fluido), res.G)[2].ravel())
        Pr.append(pr.Pr.ravel())
    Re, Pr = np.concatenate(Re), np.concatenate(Pr)
    print(f"\n--- Validez de la correlacion de Sabourishirazi " + "-" * 46)
    print(f"Re {Re.min():.4g} - {Re.max():.4g} (Tabla 4, {s.tipo}: {s.Re[0]:.0f} - "
          f"{s.Re[1]:.0f}): {np.mean((Re < s.Re[0]) | (Re > s.Re[1]))*100:.1f} % fuera")
    print(f"Pr {Pr.min():.3g} - {Pr.max():.3g} (1 - 13): "
          f"{np.mean((Pr < PR_RANGO[0]) | (Pr > PR_RANGO[1]))*100:.1f} % fuera")
    if s.fuera:
        print(f"Fuera de la Tabla 6: {', '.join(s.fuera)}  <-- EXTRAPOLACION")


def resumen(res):
    """El resumen de receptor_octogonal.py, más los datos de la correlación."""
    ro.resumen(res)
    s = res.receptor.sab
    a1, a2, a4, a5, a6, b = s.a
    print(f"\n--- Nusselt interno: Sabourishirazi et al. (2025), ec. (8) " + "-" * 35)
    print(f"Canal {s.tipo.upper()} | s_D,eff = {s.s_D_eff*1e3:.2f} mm, theta = {s.theta:.2f} grados")
    print(f"d*_ic,max = {s.d_ic:.4f}, d*_ws = {s.d_ws:.4f}, theta** = {s.theta_2:.4f}, "
          f"s_R = {s.s_R:.4f}")
    print(f"a1..a6 = {a1}, {a2}, {a4}, {a5}, {a6}, b = {b} -> Nu = {s.factor:.4g}*Re^{b}*Pr^0.4")
    validez(res)


# =============================================================================
# Caso de estudio
# =============================================================================

if __name__ == "__main__":
    mapa = MapaMitad.desde_medio_y_pico(L=ALTURA, w=ancho_placa(), n_placas=N_PLACAS // 2,
                                        q_medio=Q_MEDIO, q_pico=Q_PICO)
    principal = ReceptorSabourishirazi(mapa=mapa).resolver()
    resumen(principal)

    if BARRIDO_LAMBDA:
        resultados = [principal if lam == LAMBDA else
                      replace(principal.receptor, lam=lam).resolver(G0=principal.G)
                      for lam in sorted(set(LAMBDAS) | {LAMBDA}, reverse=True)]
        ro.tabla_lambda(resultados)

    if ro.COMPARAR_RANGO:
        ro.comparacion_rango(principal)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no esta instalado: se omiten las figuras.")
    else:
        ro.estilo_latex(plt)
        ro.dibujar_mapa(principal, plt)
        ro.dibujar_placas(principal, plt)
        ro.dibujar_correlaciones(principal, plt)
        plt.show()
