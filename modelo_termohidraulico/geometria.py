"""Geometrías de panel pillow-plate a estudiar.

Los tres paneles son los de la Tabla 1 de O. Arsenyeva et al., que a su vez los
toma del estudio experimental de M. Piper et al. Todas las longitudes en metros.

SOBRE b: la separación entre paneles NO forma parte de la Tabla 1, porque en el
método del artículo es la variable de diseño que se barre. El valor que lleva
aquí cada panel es solo un punto de partida razonable, y los modelos lo
sobrescriben caso por caso con dataclasses.replace(). Los valores actuales son
los óptimos publicados: 5.5 y 20 mm salen de la Tabla 5 (primer caso de estudio)
y 1.5 mm de la Tabla 9, variante HE1 (segundo caso de estudio).

SOBRE d_sp: tampoco está en la Tabla 1. Se toma 7.2 mm, el extremo inferior del
rango de 7.2 a 10 mm que dan Piper et al.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class PPHEGeometry:
    delta_pp: float  # Espesor de la chapa [m]
    b_i: float       # Expansión interna del panel, altura del canal interno [m]
    b: float         # Separación entre paneles, en su punto más estrecho [m]. Ver nota del módulo
    s_2l: float      # Paso longitudinal entre soldaduras, el 2*s_L de la Tabla 1 [m]
    s_t: float       # Paso transversal entre soldaduras [m]
    d_sp: float      # Diámetro del punto de soldadura [m]. Ver nota del módulo
    w_pp: float      # Anchura del panel, soldaduras de borde incluidas [m]
    l_pp: float      # Longitud del panel [m]
    w_e: float       # Anchura de la soldadura de borde [m]
    Fx: float        # Área real frente a área plana sin hidroformar [-]


# =============================================================================
# Geometrías a estudiar (Tabla 1, artículo de O. Arsenyeva et al.)
# =============================================================================

PPHE1 = PPHEGeometry(
    delta_pp=0.8e-3,
    b_i=3.4e-3,
    b=5.5e-3,
    s_2l=42e-3,
    s_t=72e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1.0,
    w_e=15e-3,
    Fx=1.010,
)

PPHE2 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=3.0e-3,
    b=1.5e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1.0,
    w_e=15e-3,
    Fx=1.007,
)

PPHE3 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=7.0e-3,
    b=20e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1.0,
    w_e=15e-3,
    Fx=1.045,
)
