"""Balance de energía del receptor, porción a porción, y temperatura de salida.

La placa L x W se parte en N filas y M columnas, y la malla la dicta el panel:
cada porción es una celda del patrón de soldaduras, dx = s_T de ancho y dy = s_L
de alto (s_L es la mitad del s_2l de la geometría). Por eso M = W/s_T y
N = L/s_L, redondeados al entero más próximo; dx y dy se reajustan luego para
cubrir la placa exacta. Con PPHE1 y la placa de 1.5 x 1.5 m salen 21 x 71.

El aire entra por abajo a T_ent uniforme, sube fila a fila y sale por y = L. La
salida de una fila es la entrada de la de encima. El resultado es el perfil
T(x) en y = L.

BALANCE DE UNA PORCIÓN. Sobre la cara exterior incide q''(x,y) del mapa
gaussiano. De ahí, una parte se pierde al ambiente por convección y otra al
cielo por radiación, ambas gobernadas por la temperatura de la pared T_w; el
resto atraviesa la chapa por conducción y pasa al fluido por convección interna:

    q_inc = h_ext*(T_w - T_amb)*A + eps*sigma*(T_w^4 - T_cielo^4)*A + q_fluido
    q_fluido = A_f*(T_w - T_f)/(e/K_ACERO + 1/h_int)      [T_f = (T_ent+T_sal)/2]
    q_fluido = G_col*(h(T_sal) - h(T_ent))                [entalpías de props_fluido]

Dos ecuaciones y dos incógnitas, T_w y T_sal, no lineales por el T_w^4 y porque
las propiedades del fluido dependen de la temperatura. Se resuelven por
bisección anidada: para cada T_sal de prueba se despeja T_w, y con él el calor
al fluido, hasta que cuadra con el salto de entalpía. La bisección se usa por
robustez: nunca diverge, y el coste no importa con mallas de este tamaño.

A y A_f NO SON LA MISMA ÁREA. La porción recibe flujo y pierde calor por toda su
cara, A = dx*dy, pero solo cede calor al fluido por A_f, lo que queda de A al
quitarle lo que no tiene aire detrás: el punto de soldadura de la celda y, en
las dos columnas de los extremos, la soldadura de borde w_e. Se supone que la
chapa está a la misma T_w en toda la porción, de modo que el calor que cae
sobre la soldadura llega al fluido por conducción lateral a través de la chapa
que la rodea; lo que se pierde es superficie de intercambio, no calor incidente.

LA SECCIÓN DE PASO es la de la ec. (16) de O. Arsenyeva, f_ch = b_i/sqrt(2) *
(W - 2*w_e), con el ancho de la placa en el papel del w_pp del panel. El gasto
se reparte entre columnas en proporción a su anchura de paso: las interiores
llevan dx y las dos de los extremos dx - w_e. Así la velocidad, u = G/(rho*f_ch)
por la ec. (3), es la misma en todas las columnas a igual temperatura.

MEZCLA ENTRE COLUMNAS, el factor LAMBDA. Cada fila se resuelve primero con las
columnas adiabáticas entre sí, lo que da una entalpía de salida h_ad,i por
columna. El caso opuesto, mezcla completa, es aquel en que toda la fila sale a
la misma temperatura, la de la entalpía media ponderada por gasto h_mez. La
salida real de la fila pondera los dos:

    h_sal,i = LAMBDA*h_ad,i + (1 - LAMBDA)*h_mez

con LAMBDA = 1 columnas adiabáticas y LAMBDA = 0 mezcla completa en cada fila.
Se pondera la entalpía y no la temperatura porque así la mezcla conserva la
energía exactamente para cualquier LAMBDA: la media ponderada de las h_sal,i es
h_mez. Con cp constante las dos ponderaciones serían la misma cosa.

PÉRDIDA DE CARGA. Es la ec. (19) del mismo artículo, integrada fila a fila
porque aquí las propiedades cambian con la temperatura, más el término de
aceleración que la (19) no lleva por tratar líquidos a propiedades constantes:

    dp = sum_j f_j*dy/de*G''^2/(2*rho_j)                  [fricción]
       + zeta_DZ/2*G''^2*(1/rho_ent + 1/rho_sal)          [zonas de distribución]
       + G''^2*(1/rho_sal - 1/rho_ent)                    [aceleración]

con G'' = G/f_ch el gasto por unidad de sección, que es el mismo en toda la
placa. El zeta_DZ de la (19) se reparte a medias entre la zona de distribución
de entrada y la de salida. Se calcula una pérdida por columna, bajando la
presión fila a fila para evaluar cada densidad a la presión local, y se da la
media ponderada por gasto.

EL CÁLCULO TÉRMICO NO DEPENDE DE LA PRESIÓN, y por eso la pérdida de carga se
calcula después y aparte. Con G'' fijo, Re = G''*de/mu no depende de rho, y cp,
k y mu son los de gas diluido (ver props_fluido): ni Nu ni h ven la presión. El
p del receptor es la presión a la ENTRADA, y solo interviene en la densidad de
los términos de la pérdida de carga.

La conductividad del acero se toma CONSTANTE (K_ACERO). La chapa es el 2.7 % de
la resistencia al fluido, de modo que su temperatura no merece una iteración
propia; ver la nota de la constante.

EL COEFICIENTE INTERNO es el del canal interno de una pillow plate, calculado
con la misma cadena que modelo_f3.py del otro paquete y con sus mismos módulos:
velocidad de la ec. (3), Re, Pr, Nu de la ec. (20) con las constantes que
asignacion deduce de la geometría, y h = Nu*k/de. La diferencia es que aquí las
propiedades del fluido no son constantes: se evalúan en cada porción a su
temperatura media con props_fluido.

HIPÓTESIS. Sin conducción lateral por la chapa entre porciones, más allá de la
que recoge LAMBDA. El reparto del gasto entre columnas es el de sus anchuras de
paso, sin redistribución por la diferencia de pérdida de carga entre ellas. No
se considera la reradiación entre porciones.
"""

from dataclasses import dataclass, field
from importlib.util import module_from_spec, spec_from_file_location
from math import pi, sqrt
from pathlib import Path

import mapa as mp
import props_fluido as pf


def _importar(nombre, ruta):
    """Carga un módulo por su ruta, sin tocar sys.path.

    Hace falta porque el otro paquete tiene un formulas.py y este también: un
    'import formulas' normal cogería el de la carpeta desde la que se ejecute,
    que es justo el error que este rodeo evita.
    """
    spec = spec_from_file_location(nombre, ruta)
    modulo = module_from_spec(spec)
    spec.loader.exec_module(modulo)
    return modulo


_TH = Path(__file__).resolve().parent.parent / "modelo_termohidraulico"
form = _importar("th_formulas", _TH / "formulas.py")
asig = _importar("th_asignacion", _TH / "asignacion.py")
geom = _importar("th_geometria", _TH / "geometria.py")


# =============================================================================
# Constantes del entorno y del material
# =============================================================================

SIGMA = 5.670374419e-8    # Constante de Stefan-Boltzmann [W/(m2*K4)]
T_AMB = 298.15            # Temperatura ambiente [K]
T_ENTORNO = 303.15        # Temperatura del entorno próximo, 30 C [K]
PESO_ENTORNO = 0.895      # Peso del entorno próximo en el sumidero radiante [-]
PESO_CIELO = 0.955        # Peso del cielo [-]

# Ecuación (6) paper M.R. Rodríguez
T_CIELO = ((PESO_ENTORNO * T_ENTORNO ** 4 + PESO_CIELO * T_AMB ** 4)
           / (PESO_ENTORNO + PESO_CIELO)) ** 0.25

H_EXT = 20.0              # Coef. de convección exterior [W/(m2*K)]. Ver nota
EPSILON = 0.85            # Emisividad de la cara expuesta [-]
ABSORTIVIDAD = 1.0        # Fracción del flujo incidente que se absorbe [-]. Ver nota
LAMBDA = 0.7              # Mezcla entre columnas: 1 adiabáticas, 0 mezcla completa [-]. Ver nota

# H_EXT: mezcla de convección natural y forzada por viento sobre una pared
# vertical caliente. La bibliografía de receptores de torre maneja de 10 W/(m2*K)
# en calma a 30 con viento; 20 es el valor intermedio habitual. Es de los
# parámetros más inciertos del modelo, pero también de los menos influyentes: a
# 1000 K la radiación se lleva del orden de cuatro veces más calor que ésta.
#
# ABSORTIVIDAD: a 1.0 el modelo es exactamente el balance pedido, con todo el
# flujo del mapa entrando en la pared. Con una pintura selectiva tipo Pyromark
# sería 0.95, y basta cambiar este número.
#
# LAMBDA: no hay dato publicado de cuánto se mezcla el flujo lateralmente en el
# canal interno de una pillow plate. Los dos extremos acotan el problema; 0.5
# es el punto medio, a falta de una simulación CFD o un ensayo que lo calibre.
# OJO: se aplica en CADA FILA, y el efecto se compone: la diferencia entre
# columnas que se crea en una fila queda multiplicada por LAMBDA en cada una de
# las siguientes. Con las 71 filas de la placa base, la dispersión del perfil de
# salida frente a la adiabática es del 72 % con LAMBDA = 0.99, el 25 % con 0.95,
# el 10 % con 0.9 y el 1 % con 0.5: el perfil sale prácticamente plano. El valor
# interesante está por encima de 0.9, y depende de s_L, que fija cuántas filas hay.


K_ACERO = 20.0            # Conductividad térmica del AISI 321 [W/(m*K)]

# K_ACERO: constante a propósito, no función de T. El acero conduce entre 16
# W/(m*K) a 100 C y 22 a 600, y la chapa del receptor se mueve, sobre todos los
# casos de trabajo ensayados, entre 400 y 1080 K de temperatura media, o sea
# entre 16 y 25 W/(m*K), con media 20.5. Se toma 20, que corresponde a unos 700
# K de chapa, el centro del rango.
#
# Que el valor exacto dé igual no es una suposición, es una cuenta: la chapa es
# el 2.7 % de la resistencia total al fluido (3.7e-5 frente a 1.3e-3 m2*K/W de
# la película interna), así que un 20 % de error en k mueve la resistencia
# total un 0.5 % y la temperatura de pared una décima de grado. Hacer k función
# de T obligaba a resolver la resistencia DENTRO de la bisección de T_w, porque
# la temperatura media de la chapa depende de la propia T_w que se busca; con k
# constante la resistencia se calcula una vez por porción y sale del bucle.
#
# El límite de servicio continuo del AISI 321 está en torno a 1150 K: si la
# pared se va por encima, el resultado deja de ser un diseño y pasa a ser un
# aviso.


def perdidas(T_w, A, h_ext=H_EXT, eps=EPSILON):
    """Calor que una porción de área A pierde al ambiente y al cielo [W].

    Está fuera de la clase porque se llama desde el interior de la bisección,
    miles de veces por porción: así se le pasan h_ext y eps como locales en vez
    de buscarlos en el objeto en cada llamada.
    """
    return (h_ext * (T_w - T_AMB) + eps * SIGMA * (T_w ** 4 - T_CIELO ** 4)) * A


def _biseccion(f, a, b, tol=1e-6, iteraciones=80):
    """Raíz de f en [a, b] por bisección. Exige que f(a) y f(b) tengan signo opuesto."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError(f"La raíz no está en [{a:.2f}, {b:.2f}]: f = {fa:.3e}, {fb:.3e}.")
    for _ in range(iteraciones):
        if b - a < tol:
            break
        m = 0.5 * (a + b)
        fm = f(m)               # Una sola evaluación por iteración: f es lo caro de aquí
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


# =============================================================================
# El receptor
# =============================================================================


@dataclass(frozen=True)
class Receptor:
    """Receptor de placa plana con el canal interno de una pillow plate."""

    mapa: mp.MapaGaussiano                  # Mapa de flujo incidente
    G: float                                # Gasto másico total [kg/s]
    T_ent: float                            # Temperatura de entrada, uniforme [K]
    p: float = 10e5                         # Presión del aire a la entrada [Pa]
    lam: float = LAMBDA                     # Mezcla entre columnas [-]
    panel: object = geom.PPHE1              # Geometría del panel pillow-plate
    fluido: object = pf.AIRE                # Fluido de trabajo
    h_ext: float = H_EXT
    eps: float = EPSILON
    alfa: float = ABSORTIVIDAD
    N: int = field(init=False)              # Porciones en altura, L/s_L
    M: int = field(init=False)              # Porciones en anchura, W/s_T
    n: object = field(init=False)           # Constantes n1..n5 de la geometría
    seccion: float = field(init=False)      # Sección de paso total del canal interno [m2]
    de: float = field(init=False)           # Diámetro equivalente del canal interno [m]
    anchos: tuple = field(init=False)       # Anchura de paso de cada columna [m]
    frac_soldadura: float = field(init=False)  # Fracción de área de puntos de soldadura [-]

    def __post_init__(self):
        """Fija de una vez lo que solo depende de la geometría, no de T."""
        i, W, L = self.panel, self.mapa.W, self.mapa.L
        fijar = lambda campo, valor: object.__setattr__(self, campo, valor)
        if not 0.0 <= self.lam <= 1.0:
            raise ValueError(f"lam tiene que estar entre 0 y 1, no {self.lam}.")

        s_L = i.s_2l / 2                                  # Paso entre filas de soldaduras
        M, N = max(round(W / i.s_t), 1), max(round(L / s_L), 1)
        dx = W / M
        if dx <= i.w_e or (M == 1 and W <= 2 * i.w_e):
            raise ValueError(f"Las soldaduras de borde ({i.w_e*1e3:.0f} mm) no caben "
                             f"en una columna de {dx*1e3:.0f} mm.")
        fijar("M", M)
        fijar("N", N)
        # Las columnas de los extremos pierden la soldadura de borde
        anchos = [dx] * M
        anchos[0] -= i.w_e
        anchos[-1] -= i.w_e
        fijar("anchos", tuple(anchos))
        # Un punto de soldadura por celda s_T x s_L: el patrón es al tresbolillo,
        # con dos puntos por cada celda s_T x s_2l
        fijar("frac_soldadura", (pi * i.d_sp ** 2 / 4) / (i.s_t * s_L))

        fijar("n", asig.asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i))
        fijar("seccion", form.f_chI(i.b_i, W, i.w_e))    # Ec. (16), w_pp = W
        fijar("de", form.deI(i.b_i))                      # Ec. (14)

    # -- Canal interno ------------------------------------------------------

    def h_interno(self, T):
        """Coeficiente de película del canal interno a la temperatura T [W/(m2*K)].

        Misma cadena que modelo_f3.py: ec. (3) -> Re, Pr -> Nu (20) -> h.
        Devuelve también u y Re, que hacen falta para comprobar que el punto de
        trabajo es razonable.
        """
        pr = self.fluido.propiedades(T, self.p)
        de, n = self.de, self.n
        u = self.G / (pr.rho * self.seccion)                 # Ec. (3)
        Re = form.Re(rho=pr.rho, u=u, dh=de, mu=pr.mu)
        Nu = form.NuI(n3=n.n3, n4=n.n4, n5=n.n5, Re=Re, Pr=pr.Pr)
        return Nu * pr.k / de, u, Re

    def friccion(self, T, p):
        """Pérdida de carga por fricción por unidad de longitud, a T y p [Pa/m].

        Primer término de la ec. (19), f/de*rho*u^2/2, con f de la ec. (18).
        """
        pr = self.fluido.propiedades(T, p)
        u = self.G / (pr.rho * self.seccion)
        Re = form.Re(rho=pr.rho, u=u, dh=self.de, mu=pr.mu)
        return form.fI(n1=self.n.n1, Re=Re, n2=self.n.n2) / self.de * pr.rho * u ** 2 / 2

    # -- Balance de una porción ----------------------------------------------

    def _porcion(self, q_inc, T_ent, G_col, A, A_f):
        """Resuelve una porción. Devuelve (T_sal, T_w, q_fluido).

        A es el área que recibe flujo y pierde calor; A_f, la que lo cede al
        fluido. Todo lo que hace falta dentro de las bisecciones se saca del
        objeto aquí, una sola vez: lo de dentro se llama miles de veces por
        porción.
        """
        fluido, e = self.fluido, self.panel.delta_pp
        h_ext, eps = self.h_ext, self.eps
        h_interno = self.h_interno
        q_abs = self.alfa * q_inc

        def cerrar(T_sal):
            """Calor al fluido y pared, para una temperatura de salida de prueba."""
            T_f = 0.5 * (T_ent + T_sal)
            h_int = h_interno(T_f)[0]
            # Conductancia chapa + película interna de toda la porción [W/K].
            # Con K_ACERO constante no depende de T_w, así que se calcula aquí,
            # una vez, en vez de en cada evaluación de la bisección de abajo.
            UA = A_f / (e / K_ACERO + 1.0 / h_int)

            def balance_pared(T_w):
                return q_abs - perdidas(T_w, A, h_ext, eps) - UA * (T_w - T_f)

            # La pared está entre el cielo (pierde más de lo que recibe) y la
            # temperatura a la que solo la radiación ya se lleva todo el flujo.
            T_w_max = max(T_f, (q_abs / (A * eps * SIGMA) + T_CIELO ** 4) ** 0.25) + 1.0
            T_w = _biseccion(balance_pared, T_CIELO, T_w_max)
            return UA * (T_w - T_f), T_w

        def desequilibrio(T_sal):
            q_fluido, _ = cerrar(T_sal)
            return q_fluido - G_col * (fluido.h(T_sal) - fluido.h(T_ent))

        # Cota superior: todo el calor absorbido al fluido, con el cp de la
        # entrada, que es el menor del tramo. Como el cp del aire crece con T, la
        # temperatura real queda por debajo.
        T_max = min(T_ent + q_abs / (G_col * fluido.cp(T_ent)), fluido.T_max)
        if desequilibrio(T_max) > 0:
            raise ValueError(
                f"El fluido se saldría de {fluido.T_max:.0f} K en una porción. "
                f"Sube el gasto o baja el flujo incidente."
            )
        T_sal = _biseccion(desequilibrio, fluido.T_min, T_max)
        q_fluido, T_w = cerrar(T_sal)
        return T_sal, T_w, q_fluido

    # -- Mezcla entre columnas -------------------------------------------------

    def _mezclar(self, T_ad, G_col):
        """Salida de una fila, ponderando adiabático y mezcla completa con lam.

        La pondera en entalpía para que la mezcla conserve la energía: ver la
        nota del módulo.
        """
        lam, h = self.lam, self.fluido.h
        T_bajo, T_alto = min(T_ad), max(T_ad)
        if lam == 1.0 or T_alto - T_bajo < 1e-9:
            return list(T_ad)
        h_ad = [h(T) for T in T_ad]
        h_mez = sum(g * hi for g, hi in zip(G_col, h_ad)) / sum(G_col)
        # Cada objetivo está entre la h mínima y la máxima de la fila; el margen
        # de 1 K es solo para que el redondeo no lo deje fuera del intervalo.
        a, b = max(T_bajo - 1.0, self.fluido.T_min), min(T_alto + 1.0, self.fluido.T_max)
        salida = []
        for h_i in h_ad:
            objetivo = lam * h_i + (1.0 - lam) * h_mez
            salida.append(_biseccion(lambda T: h(T) - objetivo, a, b, tol=1e-9))
        return salida

    # -- Recorrido de la placa -----------------------------------------------

    def resolver(self):
        """Resuelve la placa entera, fila a fila de abajo arriba."""
        M, N, porcion = self.M, self.N, self._porcion
        dy = self.mapa.L / N
        A = (self.mapa.W / M) * dy
        G_col = [self.G * a / sum(self.anchos) for a in self.anchos]
        A_f = [a * dy * (1.0 - self.frac_soldadura) for a in self.anchos]
        flujo = self.mapa.mapa_nodal(M, N)                # [W/m2] medio de cada porción

        T_fluido = []                                     # T de salida de cada fila, ya mezclada
        T_pared = []
        q_util = 0.0
        T = [self.T_ent] * M
        for j in range(N):                                # De abajo arriba
            T_ad, T_w = [0.0] * M, [0.0] * M
            for i in range(M):
                T_ad[i], T_w[i], q = porcion(flujo[j][i] * A, T[i], G_col[i], A, A_f[i])
                q_util += q
            T = self._mezclar(T_ad, G_col)
            T_fluido.append(T)
            T_pared.append(T_w)

        return Resultado(receptor=self, T_fluido=T_fluido, T_pared=T_pared, q_util=q_util,
                         G_col=G_col, **self._perdida_carga(T_fluido))

    def _perdida_carga(self, T_fluido):
        """Pérdida de carga de cada columna, por términos [Pa].

        Recorre cada columna de abajo arriba bajando la presión fila a fila,
        porque la densidad a la que se evalúa cada término es la de la presión
        local: con 100 kPa de pérdida sobre 10 bar, tomarla constante daba un
        5 % de menos. Se hace después del cálculo térmico, y aparte, porque éste
        no depende de la presión: ver la nota del módulo.
        """
        fluido, M, dy = self.fluido, self.M, self.mapa.L / self.N
        G2 = (self.G / self.seccion) ** 2                 # G''^2 [kg2/(m4*s2)]
        zeta = form.ZETA_DZ
        friccion, distribucion, aceleracion = [0.0] * M, [0.0] * M, [0.0] * M
        for i in range(M):
            T, p = self.T_ent, self.p
            v = 1.0 / fluido.rho(T, p)
            distribucion[i] = zeta / 2 * G2 * v           # Zona de distribución de entrada
            p -= distribucion[i]
            for fila in T_fluido:
                T_sal = fila[i]
                dp = self.friccion(0.5 * (T + T_sal), p) * dy
                friccion[i] += dp
                v_sal = 1.0 / fluido.rho(T_sal, p - dp)
                dp_acel = G2 * (v_sal - v)
                aceleracion[i] += dp_acel
                p -= dp + dp_acel
                if p <= 0.0:
                    raise ValueError(
                        f"La pérdida de carga agota los {self.p/1e5:.1f} bar de entrada "
                        f"antes de la salida. Sube la presión o baja el gasto."
                    )
                T, v = T_sal, v_sal
            distribucion[i] += zeta / 2 * G2 * v          # Zona de distribución de salida
        return dict(dp_friccion=friccion, dp_distribucion=distribucion,
                    dp_aceleracion=aceleracion)


@dataclass(frozen=True)
class Resultado:
    """Campo de temperaturas de la placa, pérdida de carga y balance global."""

    receptor: Receptor
    T_fluido: list          # T_fluido[j][i]: salida de la porción (columna i, fila j) [K]
    T_pared: list           # Temperatura de la pared de cada porción [K]
    q_util: float           # Calor total absorbido por el fluido [W]
    G_col: list             # Gasto de cada columna [kg/s]
    dp_friccion: list       # Pérdida de carga de cada columna, por término [Pa]
    dp_distribucion: list
    dp_aceleracion: list

    @property
    def perfil_salida(self):
        """Temperatura del fluido en y = L, columna a columna [K]. Lo pedido."""
        return self.T_fluido[-1]

    @property
    def T_salida_media(self):
        """Temperatura de mezcla del fluido a la salida, por entalpía [K]."""
        r = self.receptor
        h_media = r.fluido.h(r.T_ent) + self.q_util / r.G
        return _biseccion(lambda T: r.fluido.h(T) - h_media, r.fluido.T_min, r.fluido.T_max,
                          tol=1e-9)

    @property
    def x_centros(self):
        """Coordenada x del centro de cada columna [m]."""
        W, M = self.receptor.mapa.W, self.receptor.M
        return [W * (i + 0.5) / M for i in range(M)]

    @property
    def T_pared_max(self):
        """Temperatura máxima de la pared [K]. El límite del AISI 321 son ~1150 K."""
        return max(max(fila) for fila in self.T_pared)

    def _media(self, valores):
        """Media ponderada por el gasto de cada columna."""
        return sum(g * v for g, v in zip(self.G_col, valores)) / self.receptor.G

    @property
    def dp_columnas(self):
        """Pérdida de carga total de cada columna [Pa]."""
        return [a + b + c for a, b, c in
                zip(self.dp_friccion, self.dp_distribucion, self.dp_aceleracion)]

    @property
    def perdida_carga(self):
        """Pérdida de carga del receptor, media ponderada por gasto [Pa]."""
        return self._media(self.dp_columnas)

    @property
    def rendimiento(self):
        """Calor al fluido entre calor incidente sobre la placa [-]."""
        return self.q_util / self.receptor.mapa.potencia_total()

    def resumen(self):
        """Imprime el balance, la pérdida de carga y el perfil de salida."""
        r = self.receptor
        Q = r.mapa.potencia_total()
        T_media = self.T_salida_media
        h_int, u, Re = r.h_interno(0.5 * (r.T_ent + T_media))
        C = pf.CERO_CELSIUS
        print(f"\n{r.mapa.resumen()}")
        print(f"Aire a {r.p/1e5:.1f} bar, G = {r.G:.3f} kg/s, entrada a "
              f"{r.T_ent - C:.0f} C, lambda = {r.lam:.2f}")
        print(f"Malla {r.M} x {r.N} (dx = {r.mapa.W/r.M*1e3:.1f} mm, "
              f"dy = {r.mapa.L/r.N*1e3:.1f} mm), soldaduras {r.frac_soldadura*100:.2f} % "
              f"del area de cada porcion")
        print(f"Canal interno: seccion {r.seccion*1e4:.2f} cm2, u = {u:.1f} m/s, "
              f"Re = {Re:.3g}, h = {h_int:.0f} W/(m2*K) (a la T media del fluido)")
        print(f"\nQ incidente {Q/1e3:8.1f} kW")
        print(f"Q al fluido {self.q_util/1e3:8.1f} kW   (rendimiento "
              f"{self.rendimiento*100:.1f} %)")
        print(f"Perdidas    {(Q - self.q_util)/1e3:8.1f} kW")
        dp = self.dp_columnas
        print(f"\nPerdida de carga {self.perdida_carga/1e3:.2f} kPa: friccion "
              f"{self._media(self.dp_friccion)/1e3:.2f}, distribucion "
              f"{self._media(self.dp_distribucion)/1e3:.2f}, aceleracion "
              f"{self._media(self.dp_aceleracion)/1e3:.2f} "
              f"(entre columnas, de {min(dp)/1e3:.2f} a {max(dp)/1e3:.2f})")
        print(f"\nSalida en y = L: mezcla {T_media - C:.0f} C, "
              f"maxima {max(self.perfil_salida) - C:.0f} C, "
              f"minima {min(self.perfil_salida) - C:.0f} C")
        print(f"Pared: maxima {self.T_pared_max - C:.0f} C"
              f"{'  <-- POR ENCIMA DEL LIMITE DEL AISI 321' if self.T_pared_max > 1150 else ''}")
        print("\nPerfil de salida T(x) en y = L [C]:")
        for k in range(0, r.M, 11):                       # En tramos, que no desborde
            tramo = slice(k, k + 11)
            print("  x [m]:  " + " ".join(f"{x:6.3f}" for x in self.x_centros[tramo]))
            print("  T [C]:  " + " ".join(f"{T - C:6.0f}" for T in self.perfil_salida[tramo]))

    def dibujar(self, archivo=None, mostrar=True):
        """Campo de temperaturas del fluido y perfil de salida en y = L."""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib no esta instalado: se omite la figura.")
            return None

        r = self.receptor
        campo = [[T - pf.CERO_CELSIUS for T in fila] for fila in self.T_fluido]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
        imagen = ax1.imshow(campo, cmap=mp.CMAP, origin="lower", aspect="equal",
                            extent=(0.0, r.mapa.W, 0.0, r.mapa.L),
                            interpolation="nearest")
        ax1.set_xlabel("x, anchura [m]")
        ax1.set_ylabel("y, altura [m]")
        ax1.set_title(f"Temperatura del fluido [C], malla {r.M} x {r.N}")
        fig.colorbar(imagen, ax=ax1)

        ax2.plot(self.x_centros, [T - pf.CERO_CELSIUS for T in self.perfil_salida],
                 "o-", color="firebrick")
        ax2.set_xlim(0.0, r.mapa.W)
        ax2.set_xlabel("x, anchura [m]")
        ax2.set_ylabel("T [C]")
        ax2.set_title(f"Salida del fluido en y = L, lambda = {r.lam:.2f}")
        ax2.grid(alpha=0.3)

        fig.tight_layout()
        if archivo:
            fig.savefig(archivo, dpi=200)
            print(f"Figura guardada en {archivo}")
        if mostrar:
            plt.show()
        return fig


def tabla(resultados):
    """Tabla de consola con lo esencial de cada caso resuelto."""
    encabezado = (f"{'T entrada [C]':>14}{'dp [kPa]':>11}{'Q fluido [kW]':>15}"
                  f"{'rendimiento [%]':>17}")
    print("\n" + encabezado)
    print("-" * len(encabezado))
    for res in resultados:
        print(f"{res.receptor.T_ent - pf.CERO_CELSIUS:14.0f}{res.perdida_carga/1e3:11.2f}"
              f"{res.q_util/1e3:15.1f}{res.rendimiento*100:17.1f}")


# =============================================================================
# Caso de ejemplo
# =============================================================================

if __name__ == "__main__":
    from dataclasses import replace

    # Placa de 1.5 x 1.5 m con sigmas de un tercio del lado y 200 kW/m2 de pico:
    # 236 kW sobre la placa. El gasto es el del caso anterior de 1.0 m de ancho
    # escalado con la potencia, 0.6*1.5 kg/s, para que la pared quede en la
    # misma zona de temperaturas. 10 bar: a presión atmosférica este gasto
    # pediría del orden de 300 m/s en el canal.
    receptor = Receptor(
        mapa=mp.MapaGaussiano.desde_pico(L=1.5, W=1.5, q_pico=200e3,
                                         sigma_x=0.5, sigma_y=0.5),
        G=0.9,                           # kg/s de aire
        T_ent=pf.CERO_CELSIUS + 300,     # 300 C a la entrada, uniforme
        p=10e5,
    )

    # Barrido en la temperatura de entrada, y el detalle del caso de 300 C
    resultados = [replace(receptor, T_ent=pf.CERO_CELSIUS + T).resolver()
                  for T in (200, 300, 400, 500)]
    tabla(resultados)
    resultados[1].resumen()
    resultados[1].dibujar()
