"""Receptor final: octógono de 8 pillow plates con sal solar, resuelto por mitades.

Es el caso de estudio final del TFG. receptor.py se queda como está, con su
placa de aire, para comprobaciones y demostraciones; este módulo reutiliza de él
el balance de cada porción, las pérdidas, las constantes y la cadena de Piper y
Arsenyeva del canal interno, y cambia la geometría, el fluido y lo que se
despeja.

GEOMETRÍA. El receptor es un prisma octogonal de 8 placas pillow-plate de
ALTURA = 10.5 m. El ancho de cada placa sale de que el octógono tenga
DIAMETRO = 8.5 m. Se toma como diámetro el de la circunferencia CIRCUNSCRITA
(de vértice a vértice), que es la lectura habitual del diámetro de un polígono
regular:

    w = D*sin(pi/8) = 3.253 m        [con el inscrito, D*tan(pi/8) = 3.521 m]

Basta cambiar DIAMETRO_CIRCUNSCRITO para usar la otra lectura.

SIMETRÍA: MEDIO RECEPTOR. El mapa de flujo es simétrico respecto de un plano
vertical que pasa por dos vértices opuestos del octógono, y el fluido de cada
mitad recorre su mitad sin cruzarse con el de la otra. Se resuelve una mitad:
4 placas, un mapa desarrollado de ALTURA x 4w. En sus coordenadas, x va de 0,
el plano de simetría, a 4w, y la placa k ocupa [(k-1)*w, k*w]. Todo lo que se
da del receptor completo es el doble de lo de la mitad.

MAPA DE FLUJO. Una gaussiana sobre el receptor desarrollado (8w x ALTURA) con
el pico en el plano de simetría, a media altura. Así el pico cae en el borde de
la placa 1 de cada mitad, la del centro del receptor. Se fijan el pico,
Q_PICO = 1.2 MW/m2, y la media sobre el receptor, Q_MEDIO = 0.8 MW/m2, que por
simetría es también la de la mitad. Son dos datos para dos sigmas, así que hace
falta una hipótesis más: la campana tiene la MISMA ANCHURA RELATIVA en las dos
direcciones, sigma_x/(8w) = sigma_y/ALTURA = kappa, y kappa se despeja de
q_medio/q_pico = 2/3. Es la misma forma que usan los ejemplos de mapa.py
(sigmas proporcionales al lado). El flujo de cada porción se INTEGRA con erf,
como en mapa.py, de modo que la potencia se conserva exactamente.

RECORRIDO DEL FLUIDO. La sal entra a T_ENT = 290 C en la placa 1, la del
centro, que recibe el máximo flujo. La recorre a lo largo de su altura y sale
por el otro extremo; allí se MEZCLA COMPLETAMENTE (colector) y entra con esa
temperatura uniforme en la placa 2, y así hasta la 4. El recorrido es en
serpentín: la placa 1 sube, la 2 baja, la 3 sube y la 4 baja, que es como se
unen dos placas contiguas por un colector en un extremo. Como el mapa es
simétrico en altura, el sentido no cambia el calor absorbido; sí la pérdida de
carga hidrostática, que se da aparte.

LO QUE SE DESPEJA. La salida está fijada, T_SAL = 565 C, y la incógnita es el
GASTO G y el campo de temperaturas. Por balance global,

    G = Q_util(G) / (h(T_SAL) - h(T_ENT))

y como el calor útil depende muy poco del gasto (solo a través de la
temperatura de pared, que mueve las pérdidas), el punto fijo converge enseguida;
con la secante y el gasto_inicial(), en tres pasadas. Cada pasada resuelve las
cuatro placas.

DENTRO DE CADA PLACA, exactamente el método de receptor.py: una porción por
celda de soldadura (dx = s_T, dy = s_2l), balance de cada porción con q'' del
mapa, pérdidas por convección y radiación a T_w, chapa + película interna por
el área sin soldaduras A_f, coeficiente interno de la cadena ec. (3) -> Re, Pr
-> Nu de la ec. (20) con las constantes de Piper que da asignacion, h = Nu*k/de,
sección de paso de la ec. (16) sin soldaduras de borde y diámetro equivalente de
la (14). Las propiedades de la sal se evalúan en cada porción a su temperatura
media. La mezcla entre columnas, fila a fila, con LAMBDA:

    h_sal,i = lam*h_ad,i + (1 - lam)*h_mez

con lam = 1 columnas adiabáticas y lam = 0 mezcla completa en cada fila. Al
final de cada placa se mezcla todo, sea cual sea lam; lam cambia cómo se reparte
el calor DENTRO de las placas, y eso es lo que muestra el barrido de lambda
(BARRIDO_LAMBDA = True).

RESOLUCIÓN NUMÉRICA. Las ecuaciones son las de receptor.py, pero la malla es
unas 60 veces mayor (4 placas de 45 x 250 porciones con PPHE1), y resolver
porción a porción en Python puro llevaría minutos por cada gasto de prueba. Por
eso las M columnas de una fila se resuelven A LA VEZ con numpy: son
independientes entre sí hasta la mezcla, que se hace después. Las dos
bisecciones anidadas de receptor.py se sustituyen por Newton, que resuelve las
mismas ecuaciones en tres o cuatro evaluaciones en vez de unas veinticinco:

  - PARED. El balance q_abs - perdidas(T_w) - UA*(T_w - T_f) es decreciente y
    cóncavo en T_w (por el -T_w^4). Su derivada, -A*(h_ext + 4*eps*sigma*T_w^3)
    - UA, es estrictamente negativa para T_w > 0: no hay singularidad. Arrancando
    de T_f + q_abs/UA, donde el balance es negativo, Newton avanza hacia la raíz
    sin pasarse nunca.

  - SALIDA DE CADA PORCIÓN. El desequilibrio D(T_sal) = q_fluido - G_col*dh es
    decreciente. Se usa como derivada

        dD/dT_sal = -G_col*cp(T_sal) - 1/2 * UA*P'/(P' + UA)

    con P' = dperdidas/dT_w. El primer término es el del salto de entalpía; el
    segundo, lo que baja el calor al fluido al subir T_f con UA fijo. Los dos son
    negativos, así que la derivada no se anula nunca (cp > 1443 J/(kg*K) en la
    sal). Se desprecia solo cómo cambia UA con T_f a través de h_int, un término
    unas 1e5 veces menor que G_col*cp: la raíz es exacta (el residuo sí es el
    completo) y solo se pierde la convergencia cuadrática en la última cifra.
    Como salvaguarda se mantiene el intervalo [a, b] que encierra la raíz y, si
    un paso de Newton se sale de él, se biseca.

COMPROBACIÓN frente a la versión anterior con bisección (misma G, las cuatro
placas): las temperaturas de sal, pared y película difieren como mucho 1e-5 K,
que es la tolerancia de la bisección; el calor y la pérdida de carga de cada
placa, menos de 1e-8 en relativo. Probado con lam = 1, 0.7 y 0, con PPHE3, con
la sal a ~590 C y con un flujo del 1 % en que las pérdidas superan al flujo
absorbido y la sal se enfría. En todos, 3 evaluaciones por fila, ningún paso
de Newton fuera del intervalo, dD/dT_sal <= -3.4e3 W/K y dT_w <= -25 W/K.
El caso completo pasa de 10.4 s a 1.2 s.

PÉRDIDA DE CARGA. Ec. (19) integrada fila a fila, como en receptor.py, placa a
placa y en serie: fricción con el f de la ec. (18), zonas de distribución a la
entrada y salida de cada placa (zeta_DZ/2 cada una) y aceleración. Aparte, el
término hidrostático, rho*g*dy con signo según el sentido, que en un líquido no
es despreciable (una placa de 10.5 m son ~190 kPa) pero se compensa casi del
todo entre placas que suben y bajan.

HIPÓTESIS. Las de receptor.py (sin conducción lateral más allá de lam, gasto
igual por columna, sin reradiación entre porciones) y además: un único canal
interno por placa, que lleva todo el gasto de la mitad, y colectores ideales
entre placas (mezcla completa, sin pérdida de carga propia ni pérdidas de
calor).
"""

from dataclasses import dataclass, field, replace
from math import pi, sin, tan

import numpy as np

import mapa as mp
import props_fluido as pf
import receptor as rec

form, asig, geom = rec.form, rec.asig, rec.geom
C = pf.CERO_CELSIUS


# =============================================================================
# Datos del caso de estudio
# =============================================================================

N_PLACAS = 8                  # Placas del octógono [-]
ALTURA = 10.5                 # Altura de las placas [m]
DIAMETRO = 8.5                # Diámetro del octógono [m]
DIAMETRO_CIRCUNSCRITO = True  # True: de vértice a vértice; False: entre caras

Q_MEDIO = 0.8e6               # Flujo medio sobre el receptor [W/m2]
Q_PICO = 1.2e6                # Flujo máximo [W/m2]

T_ENT = C + 290.0             # Entrada de la sal [K]
T_SAL = C + 565.0             # Salida de la sal, impuesta [K]

LAMBDA = rec.LAMBDA           # Mezcla entre columnas del caso principal [-]
LAMBDAS = (1.0, 0.99, 0.95, 0.9, 0.7, 0.0)   # Barrido de lambda

GRAVEDAD = 9.80665            # [m/s2]
T_PELICULA_MAX = C + 600.0    # Descomposición de la sal: límite de la pared mojada [K]
T_ACERO_MAX = 1150.0          # Límite de servicio del AISI 321, el de receptor.py [K]
TOL_FILA = 1e-9               # Tolerancia de Newton en la salida de cada porción [K]

BARRIDO_LAMBDA = False        # True: resuelve también LAMBDAS y saca su tabla (x6 de tiempo)


def ancho_placa(D=DIAMETRO, n=N_PLACAS, circunscrito=DIAMETRO_CIRCUNSCRITO):
    """Lado del polígono regular de n lados y diámetro D [m]. Ver la nota del módulo."""
    return D * sin(pi / n) if circunscrito else D * tan(pi / n)


# =============================================================================
# Mapa de flujo de la mitad del receptor
# =============================================================================


@dataclass(frozen=True)
class MapaMitad:
    """Gaussiana del receptor desarrollado, vista sobre la mitad x en [0, W].

    x = 0 es el plano de simetría, donde está el pico; W = n_placas*w es la
    anchura de la mitad. La y, como en mapa.py, desde abajo. El receptor
    completo mide 2W, y la campana está centrada en (0, L/2).
    """

    L: float            # Altura [m]
    w: float            # Ancho de cada placa [m]
    n_placas: int       # Placas de la mitad [-]
    q_pico: float       # Flujo en el pico [W/m2]
    sigma_x: float      # [m]
    sigma_y: float      # [m]

    @property
    def W(self):
        """Anchura de la mitad [m]."""
        return self.n_placas * self.w

    @classmethod
    def desde_medio_y_pico(cls, L, w, n_placas, q_medio, q_pico):
        """Mapa con el pico y la media dados, sigmas proporcionales a los lados.

        sigma_x = kappa*(2W) y sigma_y = kappa*L, con kappa tal que el flujo medio
        sobre la mitad (igual al del receptor entero) sea q_medio. La relación
        media/pico crece con kappa de 0 a 1, así que la bisección basta.
        """
        if not 0.0 < q_medio < q_pico:
            raise ValueError("Hace falta 0 < q_medio < q_pico.")
        W = n_placas * w

        def exceso(kappa):
            m = cls(L, w, n_placas, q_pico, kappa * 2 * W, kappa * L)
            return m.q_medio() - q_medio

        kappa = rec._biseccion(exceso, 1e-3, 1e3, tol=1e-12)
        return cls(L, w, n_placas, q_pico, kappa * 2 * W, kappa * L)

    def _fx(self, a, b):
        return mp._fraccion_gaussiana(a, b, 0.0, self.sigma_x)

    def _fy(self, a, b):
        return mp._fraccion_gaussiana(a, b, self.L / 2.0, self.sigma_y)

    def potencia(self, x0=0.0, x1=None):
        """Potencia incidente entre x0 y x1, toda la altura [W]. Exacta, con erf."""
        x1 = self.W if x1 is None else x1
        return (self.q_pico * 2 * pi * self.sigma_x * self.sigma_y
                * self._fx(x0, x1) * self._fy(0.0, self.L))

    def q_medio(self):
        """Flujo medio sobre la mitad [W/m2]."""
        return self.potencia() / (self.W * self.L)

    def q(self, x, y):
        """Flujo en (x, y) [W/m2]. Admite arrays."""
        return self.q_pico * np.exp(-x ** 2 / (2 * self.sigma_x ** 2)
                                    - (y - self.L / 2) ** 2 / (2 * self.sigma_y ** 2))

    def mapa_nodal(self, x0, ancho, M, N):
        """Flujo MEDIO de cada porción de la franja [x0, x0 + ancho], N x M [W/m2].

        Fila j de abajo arriba y columna i creciendo en x, como mapa.mapa_nodal.
        """
        dx, dy = ancho / M, self.L / N
        fx = np.array([self._fx(x0 + i * dx, x0 + (i + 1) * dx) for i in range(M)])
        fy = np.array([self._fy(j * dy, (j + 1) * dy) for j in range(N)])
        return self.q_pico * 2 * pi * self.sigma_x * self.sigma_y * np.outer(fy, fx) / (dx * dy)

    def resumen(self):
        return (f"Mitad {self.L:.3g} x {self.W:.3g} m ({self.n_placas} placas de "
                f"{self.w:.3f} m) | sigma = ({self.sigma_x:.2f}, {self.sigma_y:.2f}) m | "
                f"q_pico = {self.q_pico/1e6:.3f} MW/m2 | q_medio = {self.q_medio()/1e6:.3f} "
                f"MW/m2 | Q = {self.potencia()/1e6:.2f} MW")


# =============================================================================
# Numérica vectorizada
# =============================================================================


def _pared(q_abs, T_f, UA, A, h_ext, eps):
    """T_w que cierra el balance de pared de cada porción, por Newton [K].

    Balance: q_abs - perdidas(T_w) - UA*(T_w - T_f) = 0, decreciente y cóncavo en
    T_w. Desde T_f + q_abs/UA el balance es negativo (ahí todo iría al fluido y
    aún quedan las pérdidas), y en una función cóncava y decreciente Newton
    avanza hacia la raíz desde ese lado sin pasarse nunca. Ver la nota del módulo.
    """
    T = T_f + q_abs / UA
    for _ in range(50):
        f = q_abs - rec.perdidas(T, A, h_ext, eps) - UA * (T - T_f)
        df = -A * (h_ext + 4 * eps * rec.SIGMA * T ** 3) - UA
        paso = f / df
        T = T - paso
        if np.max(np.abs(paso)) < 1e-8:
            return T
    raise RuntimeError("La temperatura de pared no converge.")


# =============================================================================
# Resultados
# =============================================================================


@dataclass
class Placa:
    """Resultado de una placa. Los campos van en orientación física, fila 0 abajo."""

    numero: int
    sube: bool
    x0: float                 # Borde izquierdo en el mapa de la mitad [m]
    T_ent: float              # Entrada, uniforme [K]
    T_sal: float              # Salida tras la mezcla completa del colector [K]
    perfil_salida: np.ndarray # T(x) a la salida, antes del colector [K]
    T_fluido: np.ndarray      # Salida de cada porción, ya mezclada con lam [K]
    T_pared: np.ndarray       # Cara expuesta [K]
    T_pelicula: np.ndarray    # Cara mojada, la que ve la sal [K]
    flujo: np.ndarray         # Flujo incidente medio de cada porción [W/m2]
    Q_inc: float              # [W]
    Q_util: float
    Q_conv: float
    Q_rad: float
    Q_refl: float
    dp_friccion: float        # [Pa]
    dp_distribucion: float
    dp_aceleracion: float
    dp_hidrostatica: float    # Con signo: positiva si sube

    @property
    def dp(self):
        """Pérdida de carga de la placa, sin la hidrostática [Pa]."""
        return self.dp_friccion + self.dp_distribucion + self.dp_aceleracion

    @property
    def rendimiento(self):
        return self.Q_util / self.Q_inc


@dataclass
class Resultado:
    receptor: "ReceptorOctogonal"
    G: float                  # Gasto de la mitad [kg/s]
    placas: list
    iteraciones: int

    @property
    def Q_inc(self):
        return sum(p.Q_inc for p in self.placas)

    @property
    def Q_util(self):
        return sum(p.Q_util for p in self.placas)

    @property
    def rendimiento(self):
        return self.Q_util / self.Q_inc

    @property
    def T_sal(self):
        return self.placas[-1].T_sal

    @property
    def dp(self):
        return sum(p.dp for p in self.placas)

    @property
    def dp_hidrostatica(self):
        return sum(p.dp_hidrostatica for p in self.placas)

    @property
    def T_pared_max(self):
        return max(p.T_pared.max() for p in self.placas)

    @property
    def T_pelicula_max(self):
        return max(p.T_pelicula.max() for p in self.placas)


# =============================================================================
# El receptor
# =============================================================================


@dataclass(frozen=True)
class ReceptorOctogonal:
    """Media receptor: n_placas pillow plates en serie, con sal solar."""

    mapa: MapaMitad
    T_ent: float = T_ENT
    T_sal: float = T_SAL
    lam: float = LAMBDA
    panel: object = geom.PPHE1
    fluido: object = pf.SAL_SOLAR
    h_ext: float = rec.H_EXT
    eps: float = rec.EPSILON
    alfa: float = rec.ABSORTIVIDAD
    serpentin: bool = True
    M: int = field(init=False)
    N: int = field(init=False)
    n: object = field(init=False)
    seccion: float = field(init=False)
    de: float = field(init=False)
    frac_soldadura: float = field(init=False)
    flujos: tuple = field(init=False)     # Mapa nodal de cada placa, N x M

    def __post_init__(self):
        i, w, L = self.panel, self.mapa.w, self.mapa.L
        fijar = lambda campo, valor: object.__setattr__(self, campo, valor)
        if not 0.0 <= self.lam <= 1.0:
            raise ValueError(f"lam tiene que estar entre 0 y 1, no {self.lam}.")
        M, N = max(round(w / i.s_t), 1), max(round(L / i.s_2l), 1)
        fijar("M", M)
        fijar("N", N)
        fijar("frac_soldadura", 2 * (pi * i.d_sp ** 2 / 4) / (i.s_t * i.s_2l))
        fijar("n", asig.asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i))
        fijar("seccion", form.f_chI(i.b_i, w, 0.0))       # Ec. (16), w_pp = w, w_e = 0
        fijar("de", form.deI(i.b_i))                       # Ec. (14)
        fijar("flujos", tuple(self.mapa.mapa_nodal(k * w, w, M, N)
                              for k in range(self.mapa.n_placas)))

    # -- Canal interno ------------------------------------------------------

    def h_interno(self, T, G):
        """h, u y Re del canal interno a T, con el gasto G de la placa. Admite arrays."""
        pr = self.fluido.propiedades(T)
        u = G / (pr.rho * self.seccion)                              # Ec. (3)
        Re = form.Re(rho=pr.rho, u=u, dh=self.de, mu=pr.mu)
        Nu = form.NuI(n3=self.n.n3, n4=self.n.n4, n5=self.n.n5, Re=Re, Pr=pr.Pr)
        return Nu * pr.k / self.de, u, Re

    # -- Una fila ------------------------------------------------------------

    def _fila(self, q_inc, T_ent, G_col, A, A_f, G):
        """Resuelve las M porciones de una fila a la vez.

        Mismas ecuaciones que Receptor._porcion. Devuelve T_sal (columnas
        adiabáticas), T_w, T de la cara mojada y calor al fluido, todo por
        columna.
        """
        fluido, e = self.fluido, self.panel.delta_pp
        h_ext, eps = self.h_ext, self.eps
        q_abs = self.alfa * q_inc
        h_T_ent = fluido.h(T_ent)

        def cerrar(T_sal):
            """Calor al fluido, T_w y T de película, más UA y dP/dT_w para la derivada."""
            T_f = 0.5 * (T_ent + T_sal)
            h_int = self.h_interno(T_f, G)[0]
            UA = A_f / (e / rec.K_ACERO + 1.0 / h_int)
            T_w = _pared(q_abs, T_f, UA, A, h_ext, eps)
            q_f = UA * (T_w - T_f)
            dP = A * (h_ext + 4 * eps * rec.SIGMA * T_w ** 3)
            return q_f, T_w, T_f + q_f / (A_f * h_int), UA, dP

        # Intervalo que encierra la raíz, el de receptor.py: por debajo, el
        # mínimo del fluido (D > 0); por arriba, todo el calor al fluido con el
        # cp de la entrada, el menor del tramo (D < 0). Se arranca de la cota
        # superior.
        a = np.full_like(T_ent, fluido.T_min)
        b = np.minimum(T_ent + q_abs / (G_col * fluido.cp(T_ent)), fluido.T_max)
        T_sal = b
        for it in range(30):
            q_f, T_w, T_pel, UA, dP = cerrar(T_sal)
            D = q_f - G_col * (fluido.h(T_sal) - h_T_ent)
            if it == 0 and np.any(D > 0):
                # Solo pasa si la cota superior está recortada por T_max
                raise ValueError(
                    f"La sal se saldría de {fluido.T_max - pf.CERO_CELSIUS:.0f} C en una "
                    f"porción. Sube el gasto o baja el flujo incidente.")
            # D es decreciente: si D > 0 la raíz está a la derecha de T_sal
            a = np.where(D > 0, T_sal, a)
            b = np.where(D > 0, b, T_sal)
            # Derivada: -G*cp del salto de entalpía, más lo que el calor al
            # fluido baja al subir T_f (a UA fijo). Estrictamente negativa: no
            # se anula nunca. Ver la nota del módulo.
            dD = -G_col * fluido.cp(T_sal) - 0.5 * UA * dP / (dP + UA)
            paso = -D / dD
            if np.max(np.abs(paso)) < TOL_FILA:
                return T_sal, T_w, T_pel, q_f
            nuevo = T_sal + paso
            # Salvaguarda: si Newton se sale del intervalo, se biseca
            T_sal = np.where((nuevo <= a) | (nuevo >= b), 0.5 * (a + b), nuevo)
        raise RuntimeError("La temperatura de salida de una fila no converge.")

    def _mezclar(self, T_ad):
        """Mezcla de la fila con lam, en entalpía. Columnas de igual gasto."""
        if self.lam == 1.0:
            return T_ad
        h_ad = self.fluido.h(T_ad)
        objetivo = self.lam * h_ad + (1.0 - self.lam) * h_ad.mean()
        return self.fluido.T_desde_h(objetivo)

    # -- Una placa -----------------------------------------------------------

    def _placa(self, k, T_ent, G):
        """Resuelve la placa k (0 = la central) con entrada uniforme T_ent."""
        M, N, w, L = self.M, self.N, self.mapa.w, self.mapa.L
        dy = L / N
        A = (w / M) * dy
        A_f = A * (1.0 - self.frac_soldadura)
        G_col = G / M
        sube = (k % 2 == 0) or not self.serpentin
        flujo = self.flujos[k]
        filas = range(N) if sube else range(N - 1, -1, -1)   # En el sentido del fluido

        T_fluido, T_pared, T_pel = (np.empty((N, M)) for _ in range(3))
        entradas = np.empty((N, M))                            # Entrada de cada fila
        Q_util = Q_conv = Q_rad = 0.0
        T = np.full(M, T_ent)
        for j in filas:
            entradas[j] = T
            T_ad, T_w, T_p, q_f = self._fila(flujo[j] * A, T, G_col, A, A_f, G)
            Q_util += q_f.sum()
            Q_conv += (self.h_ext * (T_w - rec.T_AMB) * A).sum()
            Q_rad += (self.eps * rec.SIGMA * (T_w ** 4 - rec.T_CIELO ** 4) * A).sum()
            T = self._mezclar(T_ad)
            T_fluido[j], T_pared[j], T_pel[j] = T, T_w, T_p

        # Colector: mezcla completa, sea cual sea lam
        T_sal = float(self.fluido.T_desde_h(self.fluido.h(T).mean()))
        Q_inc = flujo.sum() * A
        return Placa(
            numero=k + 1, sube=sube, x0=k * w, T_ent=T_ent, T_sal=T_sal,
            perfil_salida=T.copy(), T_fluido=T_fluido, T_pared=T_pared,
            T_pelicula=T_pel, flujo=flujo, Q_inc=Q_inc, Q_util=Q_util,
            Q_conv=Q_conv, Q_rad=Q_rad, Q_refl=(1.0 - self.alfa) * Q_inc,
            **self._perdida_carga(entradas, T_fluido, T, G, sube),
        )

    def _perdida_carga(self, entradas, salidas, T_ultima, G, sube):
        """Pérdida de carga de una placa, media de columnas (igual gasto) [Pa].

        Ec. (19) fila a fila como en receptor.py, más la hidrostática. entradas
        y salidas son las temperaturas de entrada y salida de cada fila.
        """
        fluido, dy = self.fluido, self.mapa.L / self.N
        G2 = (G / self.seccion) ** 2
        T_m = 0.5 * (entradas + salidas)
        pr = fluido.propiedades(T_m)
        u = G / (pr.rho * self.seccion)
        Re = form.Re(rho=pr.rho, u=u, dh=self.de, mu=pr.mu)
        f = form.fI(n1=self.n.n1, Re=Re, n2=self.n.n2)
        friccion = (f * dy / self.de * pr.rho * u ** 2 / 2).sum(axis=0)
        # Entrada uniforme; salida, la de cada columna antes del colector
        v_ent = 1.0 / fluido.rho(entradas[0 if sube else -1])
        v_sal = 1.0 / fluido.rho(T_ultima)
        distribucion = form.ZETA_DZ / 2 * G2 * (v_ent + v_sal)
        aceleracion = G2 * (v_sal - v_ent)
        hidrostatica = (1 if sube else -1) * GRAVEDAD * (pr.rho * dy).sum(axis=0)
        return dict(dp_friccion=friccion.mean(), dp_distribucion=distribucion.mean(),
                    dp_aceleracion=aceleracion.mean(), dp_hidrostatica=hidrostatica.mean())

    # -- Recorrido completo y gasto ------------------------------------------

    def recorrer(self, G):
        """Las placas en serie con el gasto G. Devuelve la lista de Placa."""
        placas, T = [], self.T_ent
        for k in range(self.mapa.n_placas):
            placas.append(self._placa(k, T, G))
            T = placas[-1].T_sal
        return placas

    def gasto_inicial(self):
        """Estimación del gasto antes de resolver, con las pérdidas aproximadas [kg/s].

        Sin pérdidas, G = alfa*Q_inc/dh, y el punto fijo arrancaba ahí, un 3.5 %
        por encima y con una pasada de más. Aquí se estiman las pérdidas placa a
        placa: la sal de cada placa a la temperatura media de su tramo
        (repartiendo el salto de entalpía según la potencia de cada placa) y la
        pared a T_f + q''*(e/k + 1/h_int), sin contar lo que las propias
        pérdidas la enfrían. Queda a pocas décimas de kg/s del gasto final.
        """
        fl, e = self.fluido, self.panel.delta_pp
        h_ent, dh = fl.h(self.T_ent), fl.h(self.T_sal) - fl.h(self.T_ent)
        G = self.alfa * self.mapa.potencia() / dh
        A = self.mapa.w * self.mapa.L / (self.M * self.N)
        Q = [f.sum() * A for f in self.flujos]
        perdidas, acumulado = 0.0, 0.0
        for Q_k, flujo in zip(Q, self.flujos):
            T_f = fl.T_desde_h(h_ent + dh * (acumulado + Q_k / 2) / sum(Q))
            acumulado += Q_k
            R = e / rec.K_ACERO + 1.0 / self.h_interno(T_f, G)[0]
            T_w = T_f + self.alfa * flujo * R / (1.0 - self.frac_soldadura)
            perdidas += rec.perdidas(T_w, A, self.h_ext, self.eps).sum()
        return (self.alfa * sum(Q) - perdidas) / dh

    def resolver(self, G0=None, tol=1e-3, iteraciones=30):
        """Gasto que lleva la sal de T_ent a T_sal, y el campo que resulta.

        Punto fijo G = Q_util(G)/dh (ver la nota del módulo), acelerado con la
        secante en cuanto hay dos pasadas: busca el cero de
        r(G) = Q_util(G)/dh - G, que es casi lineal en G. Arranca de G0 si se
        da, y si no, de gasto_inicial(). tol es sobre la temperatura de salida [K].
        """
        dh = self.fluido.h(self.T_sal) - self.fluido.h(self.T_ent)
        G = G0 or self.gasto_inicial()
        anterior = None                                       # (G, r) de la pasada previa
        for it in range(1, iteraciones + 1):
            placas = self.recorrer(G)
            if abs(placas[-1].T_sal - self.T_sal) < tol:
                return Resultado(receptor=self, G=G, placas=placas, iteraciones=it)
            r = sum(p.Q_util for p in placas) / dh - G
            if anterior is None or r == anterior[1]:
                G_nuevo = G + r                               # Punto fijo
            else:
                G_nuevo = G - r * (G - anterior[0]) / (r - anterior[1])
            anterior, G = (G, r), G_nuevo
        raise RuntimeError(f"El gasto no converge en {iteraciones} iteraciones.")


# =============================================================================
# Salida por consola
# =============================================================================


def resumen(res):
    """Parámetros de rendimiento termohidráulico del receptor."""
    r, m, i = res.receptor, res.receptor.mapa, res.receptor.panel
    G = res.G
    h_e, u_e, Re_e = r.h_interno(r.T_ent, G)
    h_s, u_s, Re_s = r.h_interno(r.T_sal, G)
    pr_e, pr_s = r.fluido.propiedades(r.T_ent), r.fluido.propiedades(r.T_sal)
    rho_m = 0.5 * (pr_e.rho + pr_s.rho)
    raya = "=" * 94

    print(f"\n{raya}\nRECEPTOR OCTOGONAL - {r.fluido.nombre.upper()}, lambda = {r.lam:.2f}\n{raya}")
    print(f"Octogono de {N_PLACAS} placas, D = {DIAMETRO} m "
          f"({'circunscrito' if DIAMETRO_CIRCUNSCRITO else 'inscrito'}): placas de "
          f"{m.L:.2f} x {m.w:.3f} m, {m.L*m.w:.2f} m2 cada una, "
          f"{N_PLACAS*m.L*m.w:.1f} m2 en total")
    print(m.resumen())
    print(f"Panel: s_T = {i.s_t*1e3:.0f} mm, s_2l = {i.s_2l*1e3:.0f} mm, b_i = "
          f"{i.b_i*1e3:.1f} mm, d_sp = {i.d_sp*1e3:.1f} mm, chapa {i.delta_pp*1e3:.1f} mm | "
          f"malla {r.M} x {r.N} por placa | soldaduras {r.frac_soldadura*100:.2f} %")
    print(f"Canal interno: seccion {r.seccion*1e4:.1f} cm2, d_e = {r.de*1e3:.2f} mm, "
          f"n1..n5 = {', '.join(f'{v:.4g}' for v in r.n)}")

    print(f"\n--- Balance termico (mitad | receptor completo) " + "-" * 44)
    print(f"Gasto de sal           {G:9.2f} kg/s | {2*G:9.2f} kg/s")
    print(f"Q incidente            {res.Q_inc/1e6:9.2f} MW   | {2*res.Q_inc/1e6:9.2f} MW")
    print(f"Q al fluido            {res.Q_util/1e6:9.2f} MW   | {2*res.Q_util/1e6:9.2f} MW")
    Q_conv = sum(p.Q_conv for p in res.placas)
    Q_rad = sum(p.Q_rad for p in res.placas)
    Q_refl = sum(p.Q_refl for p in res.placas)
    print(f"Perdidas: radiacion {Q_rad/1e6:.2f} MW, conveccion {Q_conv/1e6:.2f} MW, "
          f"reflexion {Q_refl/1e6:.2f} MW (mitad)")
    cierre = (res.Q_inc - res.Q_util - Q_conv - Q_rad - Q_refl) / res.Q_inc
    print(f"Rendimiento termico    {res.rendimiento*100:.2f} %   "
          f"(cierre del balance {cierre:+.1e})")
    print(f"Salida de la sal       {res.T_sal - C:.2f} C (objetivo {r.T_sal - C:.0f} C)")

    print(f"\n--- Canal interno " + "-" * 76)
    print(f"{'':12}{'T [C]':>8}{'rho':>9}{'mu [mPa s]':>12}{'Pr':>7}{'u [m/s]':>10}"
          f"{'Re':>10}{'h [W/m2K]':>12}")
    for nombre, T, pr, u, Re, h in (("entrada", r.T_ent, pr_e, u_e, Re_e, h_e),
                                     ("salida", r.T_sal, pr_s, u_s, Re_s, h_s)):
        print(f"{nombre:12}{T - C:8.0f}{pr.rho:9.0f}{pr.mu*1e3:12.3f}{pr.Pr:7.2f}"
              f"{u:10.2f}{Re:10.3g}{h:12.0f}")

    print(f"\n--- Por placa " + "-" * 80)
    print(f"{'placa':>5}{'sentido':>8}{'T ent':>7}{'T sal':>7}{'Q inc':>8}{'Q util':>8}"
          f"{'rend.':>7}{'q max':>7}{'Tw max':>8}{'Tw med':>8}{'Tpel max':>9}"
          f"{'dT sal':>8}{'dp':>8}")
    print(f"{'':>5}{'':>8}{'[C]':>7}{'[C]':>7}{'[MW]':>8}{'[MW]':>8}{'[%]':>7}"
          f"{'MW/m2':>7}{'[C]':>8}{'[C]':>8}{'[C]':>9}{'[K]':>8}{'[MPa]':>8}")
    for p in res.placas:
        dT = p.perfil_salida.max() - p.perfil_salida.min()
        print(f"{p.numero:5d}{'sube' if p.sube else 'baja':>8}{p.T_ent - C:7.1f}"
              f"{p.T_sal - C:7.1f}{p.Q_inc/1e6:8.2f}{p.Q_util/1e6:8.2f}"
              f"{p.rendimiento*100:7.2f}{p.flujo.max()/1e6:7.3f}{p.T_pared.max() - C:8.0f}"
              f"{p.T_pared.mean() - C:8.0f}{p.T_pelicula.max() - C:9.0f}{dT:8.2f}"
              f"{p.dp/1e6:8.3f}")
    print("  (dT sal: dispersion del perfil de salida antes del colector)")

    print(f"\n--- Limites " + "-" * 82)
    aviso = lambda x, lim: "  <-- POR ENCIMA DEL LIMITE" if x > lim else ""
    print(f"Pared expuesta maxima  {res.T_pared_max - C:6.0f} C  (AISI 321: "
          f"{T_ACERO_MAX - C:.0f} C){aviso(res.T_pared_max, T_ACERO_MAX)}")
    print(f"Pelicula de sal maxima {res.T_pelicula_max - C:6.0f} C  (descomposicion: "
          f"{T_PELICULA_MAX - C:.0f} C){aviso(res.T_pelicula_max, T_PELICULA_MAX)}")

    print(f"\n--- Perdida de carga (placas en serie) " + "-" * 55)
    fr = sum(p.dp_friccion for p in res.placas)
    di = sum(p.dp_distribucion for p in res.placas)
    ac = sum(p.dp_aceleracion for p in res.placas)
    print(f"Friccion {fr/1e6:.3f} MPa, distribucion {di/1e6:.3f} MPa, "
          f"aceleracion {ac/1e6:.4f} MPa")
    print(f"Total {res.dp/1e6:.3f} MPa = {res.dp/1e5:.1f} bar  |  hidrostatica neta "
          f"{res.dp_hidrostatica/1e3:+.2f} kPa  |  potencia de bombeo "
          f"{2*G*res.dp/rho_m/1e3:.1f} kW (receptor completo, "
          f"{2*G*res.dp/rho_m/(2*res.Q_util)*100:.2f} % del Q util)")


def tabla_lambda(resultados):
    """Cómo cambia el reparto del calor entre placas y dentro de ellas con lambda."""
    n = len(resultados[0].placas)
    print("\n" + "=" * 94)
    print("BARRIDO DE LAMBDA (mitad del receptor)")
    print("=" * 94)
    cab = (f"{'lambda':>7}{'G [kg/s]':>10}{'rend [%]':>10}"
           + "".join(f"{f'Q{k+1} [MW]':>9}" for k in range(n))
           + f"{'Tw max':>8}{'Tpel max':>9}" + f"{'dT sal P1':>10}{'dp [MPa]':>10}")
    print(cab)
    print("-" * len(cab))
    for res in resultados:
        p1 = res.placas[0]
        print(f"{res.receptor.lam:7.2f}{res.G:10.2f}{res.rendimiento*100:10.3f}"
              + "".join(f"{p.Q_util/1e6:9.3f}" for p in res.placas)
              + f"{res.T_pared_max - C:8.1f}{res.T_pelicula_max - C:9.1f}"
              + f"{p1.perfil_salida.max() - p1.perfil_salida.min():10.2f}{res.dp/1e6:10.3f}")
    print("  (Tw y Tpel en C; dT sal P1: dispersion del perfil a la salida de la placa 1, en K)")


# =============================================================================
# Figuras
# =============================================================================


def estilo_latex(plt):
    """Tipografía de las figuras a juego con la de LaTeX: Computer Modern.

    matplotlib trae la cmr10 de TeX, así que no hace falta tener LaTeX
    instalado ni usar usetex, que es mucho más lento. La cmr10 no tiene el
    signo menos de Unicode: los números de los ejes se escriben con mathtext,
    que usa la misma familia (fontset 'cm'). Los rótulos no llevan tildes
    porque la cmr10 tampoco las tiene.
    """
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
    })


def _limites_placas(ax, res, color):
    m = res.receptor.mapa
    for k in range(1, m.n_placas):
        ax.axvline(k * m.w, color=color, lw=0.8, ls="--", alpha=0.7)
    ax.set_xticks([m.w * (k + 0.5) for k in range(m.n_placas)])
    ax.set_xticklabels([f"Placa {k + 1}" for k in range(m.n_placas)])


def dibujar_mapa(res, plt):
    """Flujo incidente sobre la mitad del receptor, en inferno."""
    m = res.receptor.mapa
    x = np.linspace(0.0, m.W, 400)
    y = np.linspace(0.0, m.L, 320)
    X, Y = np.meshgrid(x, y)
    fig, ax = plt.subplots(figsize=(10, 5.2))
    im = ax.imshow(m.q(X, Y) / 1e6, cmap="inferno", origin="lower", aspect="equal",
                   extent=(0.0, m.W, 0.0, m.L), interpolation="bilinear")
    _limites_placas(ax, res, "white")
    ax.set_ylabel("$y$, altura [m]")
    ax.set_title("Flujo incidente sobre medio receptor [MW/m$^2$]  "
                 "($x = 0$: plano de simetria)")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    return fig


def dibujar_placas(res, plt):
    """Temperatura de la sal y de la pared expuesta en cada placa."""
    r, m = res.receptor, res.receptor.mapa
    n = m.n_placas
    fig, ejes = plt.subplots(2, n, figsize=(3.0 * n + 1.5, 9.5), sharey=True,
                             layout="constrained")
    campos = (("T_fluido", "Sal", mp.CMAP), ("T_pared", "Pared expuesta", mp.CMAP))
    for fila, (atributo, nombre, cmap) in enumerate(campos):
        valores = [getattr(p, atributo) - C for p in res.placas]
        vmin, vmax = min(v.min() for v in valores), max(v.max() for v in valores)
        for p, v, ax in zip(res.placas, valores, ejes[fila]):
            im = ax.imshow(v, cmap=cmap, origin="lower", aspect="auto", vmin=vmin, vmax=vmax,
                           extent=(p.x0, p.x0 + m.w, 0.0, m.L), interpolation="nearest")
            flecha = "sube" if p.sube else "baja"
            ax.set_title(f"Placa {p.numero} ({flecha})\n{nombre}: "
                         rf"{v.min():.0f} - {v.max():.0f} $^\circ$C", fontsize=9)
            if fila == len(campos) - 1:
                ax.set_xlabel("$x$ [m]")
        ejes[fila][0].set_ylabel("$y$, altura [m]")
        fig.colorbar(im, ax=list(ejes[fila]), shrink=0.9,
                     label=rf"$T$ {nombre.lower()} [$^\circ$C]")
    fig.suptitle(rf"Campos de temperatura, $\lambda$ = {r.lam:.2f}, $G$ = {res.G:.1f} kg/s "
                 f"por mitad")
    return fig


# =============================================================================
# Caso de estudio
# =============================================================================

if __name__ == "__main__":
    mapa = MapaMitad.desde_medio_y_pico(L=ALTURA, w=ancho_placa(), n_placas=N_PLACAS // 2,
                                        q_medio=Q_MEDIO, q_pico=Q_PICO)
    principal = ReceptorOctogonal(mapa=mapa).resolver()
    resumen(principal)

    if BARRIDO_LAMBDA:
        # Cada caso arranca del gasto del principal, que apenas cambia con lambda
        resultados = [principal if lam == LAMBDA else
                      replace(principal.receptor, lam=lam).resolver(G0=principal.G)
                      for lam in sorted(set(LAMBDAS) | {LAMBDA}, reverse=True)]
        tabla_lambda(resultados)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no esta instalado: se omiten las figuras.")
    else:
        estilo_latex(plt)
        dibujar_mapa(principal, plt)
        dibujar_placas(principal, plt)
        plt.show()
