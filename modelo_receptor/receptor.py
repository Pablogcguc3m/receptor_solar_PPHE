"""Balance de energía del receptor, porción a porción, y temperatura de salida.

La placa L x W se parte en N filas (dy = L/N) y M columnas (dx = W/M). Cada
columna es un canal independiente: el aire entra por abajo a T_ent uniforme,
atraviesa las N porciones de su columna y sale por y = L. La salida de una
porción es la entrada de la de encima. El resultado es el perfil T(x) en y = L.

BALANCE DE UNA PORCIÓN. Sobre la cara exterior incide q''(x,y) del mapa
gaussiano. De ahí, una parte se pierde al ambiente por convección y otra al
cielo por radiación, ambas gobernadas por la temperatura de la pared T_w; el
resto atraviesa la chapa por conducción y pasa al fluido por convección interna:

    q_inc = h_ext*(T_w - T_amb)*A + eps*sigma*(T_w^4 - T_cielo^4)*A + q_fluido
    q_fluido = A*(T_w - T_f)/(e/k_acero + 1/h_int)        [T_f = (T_ent+T_sal)/2]
    q_fluido = G_col*(h(T_sal) - h(T_ent))                [entalpías de props_fluido]

Dos ecuaciones y dos incógnitas, T_w y T_sal, no lineales por el T_w^4 y porque
las propiedades dependen de la temperatura. Se resuelven por bisección anidada:
para cada T_sal de prueba se despeja T_w, y con él el calor al fluido, hasta que
cuadra con el salto de entalpía. La bisección se usa por robustez: nunca
diverge, y el coste no importa con mallas de este tamaño.

EL COEFICIENTE INTERNO es el del canal interno de una pillow plate, calculado
con la misma cadena que modelo_f3.py del otro paquete y con sus mismos módulos:
velocidad de la ec. (3), Re, Pr, Nu de la ec. (20) con las constantes que
asignacion deduce de la geometría, y h = Nu*k/de. La diferencia es que aquí las
propiedades del fluido no son constantes: se evalúan en cada porción a su
temperatura media con props_fluido.

Nótese que la velocidad, y con ella h, no dependen de M: el gasto de una columna
y su sección de paso son ambos proporcionales a dx, así que el cociente se
cancela. Afinar la malla no cambia la física, solo la resolución del mapa.

HIPÓTESIS. Cada columna es adiabática respecto a sus vecinas: ni conducción
lateral por la chapa ni mezcla entre columnas. Es lo que hace que el perfil de
salida reproduzca la forma de la campana. Con conducción lateral el perfil
saldría más plano, de modo que éste es el caso conservador para el punto
caliente. Tampoco se considera la pérdida de carga, ni el reparto real del
caudal entre columnas (se supone uniforme), ni la reradiación entre porciones.
"""

from dataclasses import dataclass, field
from importlib.util import module_from_spec, spec_from_file_location
from math import sqrt
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
T_CIELO = 284.0           # Temperatura equivalente del cielo [K], Swinbank: 0.0552*T_amb**1.5
H_EXT = 20.0              # Coef. de convección exterior [W/(m2*K)]. Ver nota
EPSILON = 0.85            # Emisividad de la cara expuesta [-]
ABSORTIVIDAD = 1.0        # Fracción del flujo incidente que se absorbe [-]. Ver nota

# H_EXT: mezcla de convección natural y forzada por viento sobre una pared
# vertical caliente. La bibliografía de receptores de torre maneja de 10 W/(m2*K)
# en calma a 30 con viento; 20 es el valor intermedio habitual. Es de los
# parámetros más inciertos del modelo, pero también de los menos influyentes: a
# 1000 K la radiación se lleva del orden de cuatro veces más calor que ésta.
#
# ABSORTIVIDAD: a 1.0 el modelo es exactamente el balance pedido, con todo el
# flujo del mapa entrando en la pared. Con una pintura selectiva tipo Pyromark
# sería 0.95, y basta cambiar este número.


def k_acero(T):
    """Conductividad térmica del AISI 321 [W/(m*K)], T en K.

    Ajuste lineal a los valores de catálogo del acero (unos 16 W/(m*K) a 100 C y
    22 a 600 C). Los austeníticos conducen poco y además mejoran al calentarse,
    al revés que la mayoría de los metales. En el rango de trabajo la chapa es
    una resistencia pequeña frente a la convección interna, así que la linealidad
    del ajuste sobra para lo que se le pide.

    El límite de servicio continuo del AISI 321 está en torno a 1150 K: si la
    pared se va por encima, el resultado deja de ser un diseño y pasa a ser un
    aviso.
    """
    return 14.6 + 0.0127 * (T - pf.CERO_CELSIUS)


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
    p: float = 10e5                         # Presión del aire [Pa]
    N: int = 20                             # Porciones en altura
    M: int = 10                             # Porciones en anchura
    panel: object = geom.PPHE1              # Geometría del panel pillow-plate
    fluido: object = pf.AIRE                # Fluido de trabajo
    h_ext: float = H_EXT
    eps: float = EPSILON
    alfa: float = ABSORTIVIDAD
    n: object = field(init=False)           # Constantes n1..n5 de la geometría
    seccion: float = field(init=False)      # Sección de paso total del canal interno [m2]
    de: float = field(init=False)           # Diámetro equivalente del canal interno [m]

    def __post_init__(self):
        """Fija de una vez lo que solo depende de la geometría, no de T."""
        i = self.panel
        fijar = lambda campo, valor: object.__setattr__(self, campo, valor)
        fijar("n", asig.asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i))
        fijar("seccion", i.b_i / sqrt(2) * self.mapa.W)   # Ec. (16), sin soldaduras de borde
        fijar("de", form.deI(i.b_i))                      # Ec. (14)

    # -- Coeficiente de película interno -------------------------------------

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

    # -- Balance de una porción ----------------------------------------------

    def _porcion(self, q_inc, T_ent, G_col, A):
        """Resuelve una porción. Devuelve (T_sal, T_w, q_fluido).

        Todo lo que hace falta dentro de las bisecciones se saca del objeto
        aquí, una sola vez: lo de dentro se llama miles de veces por porción.
        """
        fluido, e = self.fluido, self.panel.delta_pp
        h_ext, eps = self.h_ext, self.eps
        h_interno = self.h_interno
        q_abs = self.alfa * q_inc

        def cerrar(T_sal):
            """Calor al fluido y pared, para una temperatura de salida de prueba."""
            T_f = 0.5 * (T_ent + T_sal)
            h_int = h_interno(T_f)[0]

            def resistencia(T_w):
                """Resistencia chapa + película interna, por unidad de área [m2*K/W]."""
                return e / k_acero(0.5 * (T_w + T_f)) + 1.0 / h_int

            def balance_pared(T_w):
                return (q_abs - perdidas(T_w, A, h_ext, eps)
                        - A * (T_w - T_f) / resistencia(T_w))

            # La pared está entre el cielo (pierde más de lo que recibe) y la
            # temperatura a la que solo la radiación ya se lleva todo el flujo.
            T_w_max = max(T_f, (q_abs / (A * eps * SIGMA) + T_CIELO ** 4) ** 0.25) + 1.0
            T_w = _biseccion(balance_pared, T_CIELO, T_w_max)
            return A * (T_w - T_f) / resistencia(T_w), T_w

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

    # -- Recorrido de la placa -----------------------------------------------

    def resolver(self):
        """Resuelve la placa entera, de abajo arriba y columna a columna."""
        M, N, porcion = self.M, self.N, self._porcion
        A = (self.mapa.W / M) * (self.mapa.L / N)
        G_col = self.G / M
        flujo = self.mapa.mapa_nodal(M, N)                # [W/m2] medio de cada porción

        T_fluido = [[0.0] * M for _ in range(N)]          # T de salida de cada porción
        T_pared = [[0.0] * M for _ in range(N)]
        q_util = 0.0
        for i in range(M):
            T = self.T_ent
            for j in range(N):                            # De abajo arriba
                T, T_w, q = porcion(flujo[j][i] * A, T, G_col, A)
                T_fluido[j][i], T_pared[j][i] = T, T_w
                q_util += q
        return Resultado(receptor=self, T_fluido=T_fluido, T_pared=T_pared, q_util=q_util)


@dataclass(frozen=True)
class Resultado:
    """Campo de temperaturas de la placa y balance global."""

    receptor: Receptor
    T_fluido: list      # T_fluido[j][i]: salida de la porción (columna i, fila j) [K]
    T_pared: list       # Temperatura de la pared de cada porción [K]
    q_util: float       # Calor total absorbido por el fluido [W]

    @property
    def perfil_salida(self):
        """Temperatura del fluido en y = L, columna a columna [K]. Lo pedido."""
        return self.T_fluido[-1]

    @property
    def x_centros(self):
        """Coordenada x del centro de cada columna [m]."""
        W, M = self.receptor.mapa.W, self.receptor.M
        return [W * (i + 0.5) / M for i in range(M)]

    @property
    def T_pared_max(self):
        """Temperatura máxima de la pared [K]. El límite del AISI 321 son ~1150 K."""
        return max(max(fila) for fila in self.T_pared)

    @property
    def rendimiento(self):
        """Calor al fluido entre calor incidente sobre la placa [-]."""
        return self.q_util / self.receptor.mapa.potencia_total()

    def resumen(self):
        """Imprime el balance y el perfil de salida."""
        r = self.receptor
        Q = r.mapa.potencia_total()
        h_int, u, Re = r.h_interno(0.5 * (r.T_ent + sum(self.perfil_salida) / r.M))
        T_media = sum(self.perfil_salida) / r.M
        print(f"\n{r.mapa.resumen()}")
        print(f"Aire a {r.p/1e5:.1f} bar, G = {r.G:.3f} kg/s, entrada a "
              f"{r.T_ent - pf.CERO_CELSIUS:.0f} C, malla {r.M} x {r.N}")
        print(f"Canal interno: u = {u:.1f} m/s, Re = {Re:.3g}, h = {h_int:.0f} W/(m2*K) "
              f"(a la T media del fluido)")
        print(f"\nQ incidente {Q/1e3:8.1f} kW")
        print(f"Q al fluido {self.q_util/1e3:8.1f} kW   (rendimiento "
              f"{self.rendimiento*100:.1f} %)")
        print(f"Perdidas    {(Q - self.q_util)/1e3:8.1f} kW")
        print(f"\nSalida en y = L: media {T_media - pf.CERO_CELSIUS:.0f} C, "
              f"maxima {max(self.perfil_salida) - pf.CERO_CELSIUS:.0f} C, "
              f"minima {min(self.perfil_salida) - pf.CERO_CELSIUS:.0f} C")
        print(f"Pared: maxima {self.T_pared_max - pf.CERO_CELSIUS:.0f} C"
              f"{'  <-- POR ENCIMA DEL LIMITE DEL AISI 321' if self.T_pared_max > 1150 else ''}")
        print("\nPerfil de salida T(x) en y = L [C]:")
        print("  x [m]:  " + " ".join(f"{x:6.3f}" for x in self.x_centros))
        print("  T [C]:  " + " ".join(f"{T - pf.CERO_CELSIUS:6.0f}" for T in self.perfil_salida))

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
        ax2.set_title("Salida del fluido en y = L")
        ax2.grid(alpha=0.3)

        fig.tight_layout()
        if archivo:
            fig.savefig(archivo, dpi=200)
            print(f"Figura guardada en {archivo}")
        if mostrar:
            plt.show()
        return fig


# =============================================================================
# Caso de ejemplo
# =============================================================================

if __name__ == "__main__":
    # Punto de trabajo escogido para que la pared no rebase el límite del AISI
    # 321 y el Reynolds quede cerca del de los casos con los que se validó la
    # correlación de Nusselt. Con 300 kW/m2 de pico y este mismo gasto, la pared
    # se va a 946 C y el acero no aguanta: hay que subir el gasto a 0.6 kg/s.
    receptor = Receptor(
        mapa=mp.MapaGaussiano.desde_pico(L=1.5, W=1.0, q_pico=200e3,
                                         sigma_x=1 / 3, sigma_y=0.5),
        G=0.35,                          # kg/s de aire
        T_ent=pf.CERO_CELSIUS + 300,     # 300 C a la entrada, uniforme
        p=10e5,                          # 10 bar: a presión atmosférica este
                                         # gasto pediría 300 m/s en el canal
        N=20, M=10,
    )
    resultado = receptor.resolver()
    resultado.resumen()
    resultado.dibujar()
