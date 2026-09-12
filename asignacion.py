"""Este script asigna los valores para las leyes de potencia de factor de fricción y Nu según el artículo de M.Piper, tablas 2 y 3"""

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

