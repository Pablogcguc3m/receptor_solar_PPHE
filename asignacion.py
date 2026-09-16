"""Constantes de las leyes de potencia del canal interno, deducidas de la geometría.

El canal interno de una pillow plate tiene un factor de fricción y un número de
Nusselt que se ajustan a leyes de potencia,

    zeta = n1 * Re**n2            (ec. 18 de O. Arsenyeva)
    Nu   = n3 * Re**n4 * Pr**n5   (ec. 20)

y cuyos coeficientes dependen de la geometría del panel. Este módulo los obtiene
de las tablas 2 y 3 del artículo de M. Piper, que las dan en función de tres
relaciones adimensionales entre los pasos de soldadura y la altura del canal.

Las tres familias cubiertas corresponden a paneles de celda longitudinal,
cuadrada y transversal. Fuera de esos rangos no hay correlación publicada, y la
función avisa en lugar de extrapolar.

NOTA: para las tres geometrías concretas del caso de estudio (PPHE1, PPHE2 y
PPHE3), el artículo de O. Arsenyeva tabula directamente sus propios coeficientes
en sus tablas 2 y 3, y NO coinciden con los que salen de aquí: hay diferencias de
hasta el 51 % en zeta y el 24 % en Nu, sobre todo en PPHE3. Para reproducir sus
resultados hay que usar los suyos; para una geometría nueva, como la del
receptor solar, esta vía es la única disponible.
"""

from typing import NamedTuple


class ConstantesPotencia(NamedTuple):
    """Coeficientes de las dos leyes de potencia del canal interno.

    Se puede desempaquetar e indexar como una lista, de modo que n[0] es n1.
    """
    n1: float   # Factor de fricción: zeta = n1*Re**n2
    n2: float
    n3: float   # Nusselt: Nu = n3*Re**n4*Pr**n5
    n4: float
    n5: float


# Rangos de validez de cada familia, tal y como los dan las tablas de M. Piper:
# (nombre, (razon_pasos min, max), (razon_soldadura min, max), (razon_altura min, max))
_FAMILIAS = (
    ("celda longitudinal", (0.57, 0.59), (0.10, 0.14), (0.042, 0.083)),
    ("celda cuadrada",     (0.99, 1.01), (0.17, 0.24), (0.071, 0.143)),
    ("celda transversal",  (1.70, 1.72), (0.17, 0.24), (0.071, 0.170)),
)


def asignacion_de_constantes(sT, s2L, dsp, h):
    """Devuelve las constantes n1..n5 que corresponden a una geometría de panel.

    Argumentos, todos en metros:
        sT   paso transversal entre soldaduras
        s2L  paso longitudinal entre soldaduras (el 2*s_L de las tablas)
        dsp  diámetro del punto de soldadura
        h    expansión interna del panel, b_i

    Lanza ValueError si la geometría cae fuera de los rangos ajustados.
    """
    razon_pasos = s2L / sT
    razon_soldadura = dsp / sT
    razon_altura = h / sT

    # Abreviaturas con los nombres de las tablas originales, para que las
    # expresiones de abajo se puedan cotejar de un vistazo con el artículo.
    a, b, c = razon_pasos, razon_soldadura, razon_altura

    if 0.57 <= a <= 0.59 and 0.10 <= b <= 0.14 and 0.042 <= c <= 0.083:
        # Celda longitudinal
        n1 = 8.74 * b + (17 * c + 0.73)
        n2 = -0.38
        n3 = 0.0775 * b + (0.38 * c + 0.005)
        n4 = 0.75
        n5 = 0.4
    elif 0.99 <= a <= 1.01 and 0.17 <= b <= 0.24 and 0.071 <= c <= 0.143:
        # Celda cuadrada
        n1 = -15.3 * b + (1.4 * c + 5.4)
        n2 = 1.725 * b + (1.11 * c - 0.66)
        n3 = 0.03 * b + (0.76 * c - 0.032)
        n4 = -1.12 * c + 0.905
        n5 = 0.4
    elif 1.70 <= a <= 1.72 and 0.17 <= b <= 0.24 and 0.071 <= c <= 0.170:
        # Celda transversal
        n1 = 1.35 * b + (2.8 * c + 0.92)
        n2 = 0.3 * b + (0.53 * c - 0.29)
        n3 = -0.163 * b + (0.711 * c + 0.022)
        n4 = 0.29 * b + (-c + 0.8)
        n5 = 0.4
    else:
        rangos = "\n".join(
            f"    {nombre:<20} s2L/sT en [{ra[0]}, {ra[1]}], "
            f"dsp/sT en [{rb[0]}, {rb[1]}], b_i/sT en [{rc[0]}, {rc[1]}]"
            for nombre, ra, rb, rc in _FAMILIAS
        )
        raise ValueError(
            f"La geometría queda fuera de las correlaciones publicadas:\n"
            f"    s2L/sT = {a:.3f}, dsp/sT = {b:.3f}, b_i/sT = {c:.3f}\n"
            f"Rangos admitidos:\n{rangos}"
        )

    return ConstantesPotencia(n1, n2, n3, n4, n5)
