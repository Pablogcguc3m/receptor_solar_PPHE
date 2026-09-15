"""Este script recoge las geometrías de los PPHE a estudiar."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PPHEGeometry:
    delta_pp: float  # Espesor de placa [m]
    b_i: float       # Altura/profundidad interna [m]
    b: float         # Distancia entre paneles, en su punto más estrecho [m].
                     # OJO: no está en la Tabla 1 del artículo; es su variable de
                     # diseño. Estos valores son los óptimos que da su Tabla 5.
    s_2l: float      # Paso longitudinal [m]
    s_t: float       # Paso transversal [m]
    d_sp: float      # Diagonal de la soldadura [m]
    w_pp: float      # Anchura de la placa [m]
    l_pp: float      # Longitud de la placa [m]
    w_e: float       # Anchura de soldadura de borde [m]
    Fx: float        # Ratio de área real con respecto al área plana inicial sin hidroformado [-]


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
    l_pp=1000e-3,
    w_e=15e-3,
    Fx=1.01
)
PPHE2 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=3.0e-3,
    b=7.5e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3,
    Fx=1.007
)
PPHE3 = PPHEGeometry(
    delta_pp=1.0e-3,
    b_i=7.0e-3,
    b=20e-3,
    s_2l=72e-3,
    s_t=42e-3,
    d_sp=7.2e-3,
    w_pp=300e-3,
    l_pp=1000e-3,
    w_e=15e-3,
    Fx=1.045
)