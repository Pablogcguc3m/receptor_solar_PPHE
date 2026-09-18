"""Fórmulas del modelo termohidráulico de intercambiadores pillow-plate (PPHE).

Salvo donde se indique otra cosa, las correlaciones proceden de

    O. Arsenyeva, J. Tran, M. Piper, E. Kenig, "An approach for pillow plate heat
    exchangers design for single-phase applications", Appl. Therm. Eng. 147
    (2019) 579-591.

El número entre paréntesis al final de cada docstring es el de la ecuación en ese
artículo.

CONVENIO DE SUBÍNDICES. En este módulo, I es el canal INTERNO (el que queda dentro
del panel soldado) y E el EXTERNO (el que forman dos paneles contiguos). Donde
aparecen 1 y 2 son los subíndices del propio artículo, en los que 1 es la
corriente CALIENTE y 2 la FRÍA. Cuidado al mezclarlos: los dos criterios solo
coinciden cuando la corriente caliente circula por el canal interno.

DESVIACIONES RESPECTO AL ARTÍCULO. Hay dos, ambas documentadas en el docstring de
la función correspondiente y justificadas en los informes de informes/:
  - wI usa 8 donde el artículo imprime 8*sqrt(2), ec. (28).
  - NTU toma valores absolutos, para admitir que el fluido 1 sea el frío.
"""

from math import log, sqrt

# =============================================================================
# Constantes del modelo
# =============================================================================

LAMBDA_W = 16.0    # Conductividad térmica de la pared [W/(m*K)]
ZETA_DZ = 1.5      # Resistencia hidráulica local de las zonas de distribución [-]
PSI_E = 0.58       # Coeficiente de corrección del canal externo, ec. (21) [-]

# Ley de potencia del factor de fricción del canal EXTERNO, ec. (18) con la
# Tabla 2. El artículo solo la ajustó para la geometría PPHE3, y la aplica a
# todas: puede asumirse independiente de la geometría del panel.
A_FRICCION_E = 2.187
N_FRICCION_E = 0.356


# =============================================================================
# Números adimensionales
# =============================================================================

def Re(rho, u, dh, mu):
    """Número de Reynolds."""
    return (rho * u * dh) / mu


def Pr(cp, mu, k):
    """Número de Prandtl."""
    return (cp * mu) / k


def NTU(t1in, t1out, t2in, t2out):
    """Número de unidades de transferencia MÍNIMO, referido al fluido 1. (23)

    Contracorriente. Se toman valores absolutos para que sea indiferente que el
    fluido 1 sea el caliente o el frío: los saltos terminales son |t1in-t2out|,
    en el extremo por el que entra el fluido 1, y |t1out-t2in| en el otro. Tal y
    como está impresa en el artículo, la expresión solo admite que el fluido 1
    sea el caliente, y con el frío el logaritmo recibe un argumento negativo.
    """
    dT_a = abs(t1in - t2out)
    dT_b = abs(t1out - t2in)
    if abs(dT_a - dT_b) < 1e-9 * max(dT_a, dT_b):
        LMTD = (dT_a + dT_b)/2      # Límite dT_a -> dT_b, donde log(1)=0 indetermina la expresión general
    else:
        LMTD = (dT_a - dT_b)/log(dT_a/dT_b)
    return abs(t1in - t1out)/LMTD


# =============================================================================
# Factores de fricción y números de Nusselt
# =============================================================================

def fI(n1, Re, n2):
    """Factor de fricción del canal INTERNO, ley de potencia. (18)

    n1 y n2 salen de asignacion.asignacion_de_constantes, que los deduce de la
    geometría del panel según las tablas de M. Piper.
    """
    return n1 * Re ** n2


def fE(Re, A=A_FRICCION_E, n=N_FRICCION_E):
    """Factor de fricción del canal EXTERNO, ley de potencia. (18)"""
    return A * Re ** (-n)


def NuI(n3, n4, n5, Re, Pr):
    """Número de Nusselt del canal INTERNO, ley de potencia. (20)

    Es la ecuación (6) del artículo de M. Piper. n3, n4 y n5 salen de
    asignacion.asignacion_de_constantes.
    """
    return n3 * Re ** n4 * Pr ** n5


def NuE(Re, Pr, f, psi=PSI_E):
    """Número de Nusselt del canal EXTERNO, correlación experimental. (21)

    f es el factor de fricción del canal externo, el que devuelve fE.
    """
    return (psi*f/8*Re*Pr)/(1.07 + 12.7*(psi*f/8)**(1/2)*(Pr**(2/3) - 1))


# =============================================================================
# Geometría del canal: diámetros equivalentes y secciones de paso
# =============================================================================

def deI(b_i):
    """Diámetro equivalente del canal INTERNO. (14)

    Es dos veces el hueco medio del canal, que vale b_i/sqrt(2).
    """
    return 2 * (b_i / sqrt(2))


def deE(b_i, b):
    """Diámetro equivalente del canal EXTERNO. (15)

    Es dos veces el hueco medio, que es el paso total (b_i+b) menos el espesor
    medio que ocupa el panel, b_i/sqrt(2).
    """
    return 2 * ((b_i + b) - b_i/sqrt(2))


def f_chI(b_i, w_pp, w_e):
    """Sección de paso de UN canal INTERNO. (16)

    Es deI/2, el hueco medio, por la anchura útil w_pp-2*w_e, que descuenta las
    soldaduras de borde.
    """
    return (b_i/sqrt(2)) * (w_pp - 2*w_e)


def f_chE(b_i, b, w_pp, w_e):
    """Sección de paso de UN canal EXTERNO. (17) Análoga a f_chI, con deE/2."""
    return ((b_i + b) - b_i/sqrt(2)) * (w_pp - 2*w_e)


# =============================================================================
# Velocidades
# =============================================================================

def velocidad_ec3(G, rho, N, f_ch):
    """Velocidad media a partir del gasto másico. (3)

    G es el gasto de toda la corriente y N el número de canales entre los que se
    reparte, de modo que N*f_ch es la sección de paso total de ese lado.

    Es la vía directa: impone el caudal y deja la pérdida de carga como
    resultado. La alternativa es wI, que hace justo lo contrario.
    """
    return G/(rho * N * f_ch)


def wI(del_P, rho, f, NTU, cp, w, U, Fx, zeta_DZ=ZETA_DZ):
    """Velocidad del canal INTERNO agotando la pérdida de carga admisible. (28)

    Ecuación implícita en w: hay que iterarla (ver iteracion_vel.py).

    OJO, DESVIACIÓN RESPECTO AL ARTÍCULO: aquí el denominador es 8 donde el
    artículo imprime 8*sqrt(2). Esa constante es incompatible con el resto de sus
    propias ecuaciones: igualando las dos expresiones de L_F (LF_ploss, de la ec.
    22, y LF_thermal, de la 26) y sustituyendo deI = 2*b_i/sqrt(2), el b_i se
    cancela y el denominador queda 8. Con 8*sqrt(2) las dos longitudes difieren
    siempre en un factor sqrt(2); con 8 coinciden exactamente.
    """
    return sqrt(del_P/rho * 1/(zeta_DZ + (f*NTU*cp*w*rho)/(8*U*Fx)))


def wE(w_I, cp1, cp2, rho1, rho2, t1in, t1out, t2in, t2out, de1, de2, wpp, we):
    """Velocidad del otro canal, por balance de masa y energía. (30)

    Los subíndices 1 y 2 son los de los argumentos que se pasen, no los del
    artículo: 1 es el canal cuya velocidad se conoce (w_I) y 2 aquel cuya
    velocidad se busca. El cociente de secciones va con la del canal 1 en el
    numerador; invertirlo rompe la conservación de la masa.
    """
    fch1 = de1/2 * (wpp - 2*we)
    fch2 = de2/2 * (wpp - 2*we)
    return w_I * (cp1*rho1*(t1in - t1out)*fch1)/(cp2*rho2*(t2out - t2in)*fch2)


# =============================================================================
# Coeficiente global, longitudes y pérdida de carga
# =============================================================================

def U(h1, h2, delta_w, lambda_w=LAMBDA_W):
    """Coeficiente global de transmisión de calor. (29)

    h1 y h2 son los coeficientes de película de los dos canales, y delta_w el
    espesor de la chapa.
    """
    return 1/(1/h1 + 1/h2 + delta_w/lambda_w)


def LF_thermal(b_i, NTU, cp, w, rho, U, Fx):
    """Longitud de canal necesaria para transferir el calor requerido. (26)

    w, cp y rho son los del fluido del canal INTERNO.
    """
    return (b_i*NTU*cp*w*rho)/(2*sqrt(2)*U*Fx)


def LF_ploss(de1, f, del_P, rho, w, zeta_DZ=ZETA_DZ):
    """Longitud de canal que consume una pérdida de carga dada. (22)

    Es la inversa de perdida_de_carga.
    """
    return 2*de1/f * (del_P/(rho*w**2) - zeta_DZ)


def perdida_de_carga(f, L_F, d_e, rho, w, zeta_DZ=ZETA_DZ):
    """Pérdida de carga del canal, zonas de distribución incluidas. (19)

    Es la inversa de LF_ploss: aquélla despeja la longitud conocida la pérdida de
    carga, y ésta la pérdida de carga conocida la longitud.
    """
    return f*L_F/d_e * rho*w**2/2 + zeta_DZ*rho*w**2
