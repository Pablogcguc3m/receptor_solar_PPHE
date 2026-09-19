"""Mapa de flujo solar incidente sobre la placa del receptor: q''(x, y) gaussiano.

Genera la función f(x, y) que da el calor incidente por unidad de superficie
sobre una placa plana de altura L y anchura W, repartido como una campana de
Gauss centrada en el centro de la placa:

    q''(x, y) = q_pico * exp( -( (x-xc)^2/(2*sx^2) + (y-yc)^2/(2*sy^2) ) )

con (xc, yc) = (W/2, L/2). Es la forma que se le supone a la mancha solar de un
campo de helióstatos sobre un receptor plano, y la que usa la práctica totalidad
de la bibliografía de receptores cuando no se dispone de un trazado de rayos.

SISTEMA DE COORDENADAS. Origen en la esquina INFERIOR IZQUIERDA de la placa
vista de frente, x hacia la derecha y creciendo con la anchura, y hacia arriba y
creciendo con la altura:

    x en [0, W]   anchura  [m]
    y en [0, L]   altura   [m]

Así el eje y es el de la dirección en que corre el fluido en un panel vertical,
que es la que luego se discretiza en el modelo nodal.

UNIDADES. SI: longitudes en m, flujo q'' en W/m2, potencias en W.

DOS FORMAS DE FIJAR LA CAMPANA. O se da el flujo de pico, que suele ser el
límite de diseño del material del receptor, o se da la potencia total que el
campo deja sobre la placa, que es lo que sale del dimensionado del campo:

    MapaGaussiano.desde_pico(L, W, q_pico, sigma_x, sigma_y)
    MapaGaussiano.desde_potencia(L, W, Q, sigma_x, sigma_y)

CUIDADO CON LA POTENCIA: LA PLACA ES FINITA. La integral de la gaussiana en todo
el plano es 2*pi*sx*sy*q_pico, pero la placa solo recoge la parte central de la
campana; el resto se pierde por los bordes (derrame o spillage). Por eso aquí la
potencia se integra SOBRE LA PLACA, con la función error:

    Q = 2*pi*sx*sy*q_pico * erf(W/(2*sqrt(2)*sx)) * erf(L/(2*sqrt(2)*sy))

Los dos erf son exactamente la fracción de campana interceptada en cada
dirección, y su producto es el rendimiento de interceptación del mapa; lo
devuelve factor_interceptacion(). Tomar 2*pi*sx*sy*q_pico como potencia de la
placa es el error clásico: con sigmas de un tercio del lado se pasa de largo un
33 %, y con sigmas del orden del lado, un 357 %.

REPARTO POR NODOS. mapa_nodal() no muestrea la gaussiana en el centro de cada
celda: integra la campana dentro de cada celda, también con erf, y devuelve el
flujo MEDIO de cada una. Con eso la suma de las potencias nodales reproduce la
potencia total exactamente, cosa que el muestreo puntual no hace: con mallas
bastas se deja por el camino varios puntos porcentuales, y el balance de energía
del receptor no cierra.

LO QUE ESTE MÓDULO NO ES. Es una distribución geométrica de flujo INCIDENTE: no
sabe de absortividad, de pérdidas por reradiación ni de convección al ambiente.
Todo eso va en el balance del receptor, aguas abajo, y toma este mapa como dato.
"""

from dataclasses import dataclass
from math import erf, exp, log, pi, sqrt

# =============================================================================
# Constantes
# =============================================================================

# Un pico de FWHM sigmas: relación entre anchura a media altura y desviación
# típica de una gaussiana, FWHM = 2*sqrt(2*ln2)*sigma.
FWHM_POR_SIGMA = 2.0 * sqrt(2.0 * log(2.0))

# Escala de color de las figuras: 'hot' recorre negro -> rojo -> amarillo ->
# blanco, es decir, el mínimo en negro y el máximo en blanco, que es como se ve
# un cuerpo incandescente y como se presentan los mapas de flujo en la
# bibliografía de receptores. 'afmhot' es la misma escala más clara y 'inferno'
# la versión perceptualmente uniforme, que llega a amarillo pero no a blanco.
CMAP = "hot"


def sigma_desde_fwhm(fwhm):
    """Convierte la anchura a media altura de la mancha en desviación típica [m].

    La bibliografía y las medidas de campo suelen dar la mancha por su FWHM o
    por su diámetro, no por sigma.
    """
    return fwhm / FWHM_POR_SIGMA


def fwhm_desde_sigma(sigma):
    """La conversión inversa de sigma_desde_fwhm [m]."""
    return FWHM_POR_SIGMA * sigma


def _fraccion_gaussiana(a, b, centro, sigma):
    """Integral de la gaussiana normalizada entre a y b, dividida por sigma*sqrt(2*pi).

    Es la probabilidad acumulada entre a y b, o sea, la fracción de campana que
    cae en ese tramo. Escrita con erf para que sea exacta y no haya que integrar
    numéricamente en ningún sitio del módulo.
    """
    t = sigma * sqrt(2.0)
    return 0.5 * (erf((b - centro) / t) - erf((a - centro) / t))


# =============================================================================
# El mapa
# =============================================================================


@dataclass(frozen=True)
class MapaGaussiano:
    """Distribución gaussiana de flujo incidente sobre una placa L x W.

    Se construye con desde_pico o desde_potencia, no directamente: así queda
    explícito cuál de los dos datos es el que se impone.
    """

    L: float          # Altura de la placa [m]
    W: float          # Anchura de la placa [m]
    q_pico: float     # Flujo en el centro de la placa [W/m2]
    sigma_x: float    # Desviación típica de la campana en anchura [m]
    sigma_y: float    # Desviación típica de la campana en altura [m]

    def __post_init__(self):
        for nombre in ("L", "W", "sigma_x", "sigma_y"):
            if getattr(self, nombre) <= 0:
                raise ValueError(f"{nombre} tiene que ser positivo.")
        if self.q_pico < 0:
            raise ValueError("q_pico no puede ser negativo.")

    # -- Constructores -------------------------------------------------------

    @classmethod
    def desde_pico(cls, L, W, q_pico, sigma_x, sigma_y=None):
        """Mapa a partir del flujo de pico [W/m2]. sigma_y = sigma_x si se omite."""
        return cls(L=L, W=W, q_pico=q_pico, sigma_x=sigma_x,
                   sigma_y=sigma_x if sigma_y is None else sigma_y)

    @classmethod
    def desde_potencia(cls, L, W, Q, sigma_x, sigma_y=None):
        """Mapa a partir de la potencia total SOBRE LA PLACA, Q [W].

        Despeja q_pico de la integral sobre la placa, no de la integral en todo
        el plano: el q_pico que sale es por tanto mayor que Q/(2*pi*sx*sy), y
        tanto más cuanto más ancha sea la campana frente a la placa.
        """
        sy = sigma_x if sigma_y is None else sigma_y
        provisional = cls(L=L, W=W, q_pico=1.0, sigma_x=sigma_x, sigma_y=sy)
        return cls(L=L, W=W, q_pico=Q / provisional.potencia_total(),
                   sigma_x=sigma_x, sigma_y=sy)

    # -- Geometría -----------------------------------------------------------

    @property
    def centro(self):
        """Centro de la placa, donde se sitúa el pico de la campana (xc, yc) [m]."""
        return (self.W / 2.0, self.L / 2.0)

    @property
    def area(self):
        """Área de la placa, L*W [m2]."""
        return self.L * self.W

    def _verificar(self, x, y):
        """Comprueba que el punto cae dentro de la placa.

        Se prefiere el error a devolver la cola de la gaussiana fuera de la
        placa: un punto fuera casi siempre es un fallo de indexado de la malla,
        y el derrame ya se calcula aparte y de forma exacta.
        """
        if not (0.0 <= x <= self.W and 0.0 <= y <= self.L):
            raise ValueError(
                f"El punto ({x:.4g}, {y:.4g}) m cae fuera de la placa "
                f"[0, {self.W:.4g}] x [0, {self.L:.4g}] m."
            )

    # -- La función f(x, y) --------------------------------------------------

    def q(self, x, y):
        """Flujo incidente en el punto (x, y) de la placa [W/m2].

        Ésta es la f(x, y) del mapa: x es la coordenada en la anchura y la de la
        altura, ambas medidas desde la esquina inferior izquierda.
        """
        self._verificar(x, y)
        xc, yc = self.centro
        return self.q_pico * exp(
            -((x - xc) ** 2 / (2.0 * self.sigma_x ** 2)
              + (y - yc) ** 2 / (2.0 * self.sigma_y ** 2))
        )

    def __call__(self, x, y):
        """Permite usar el mapa como si fuera la propia función: mapa(x, y)."""
        return self.q(x, y)

    # -- Integrales sobre la placa -------------------------------------------

    def potencia_total(self):
        """Potencia incidente sobre la placa [W]. Integral exacta, con erf.

        Es la integral de la campana DENTRO de la placa, no en todo el plano:
        ver la nota del módulo.
        """
        xc, yc = self.centro
        fx = _fraccion_gaussiana(0.0, self.W, xc, self.sigma_x)
        fy = _fraccion_gaussiana(0.0, self.L, yc, self.sigma_y)
        return self.q_pico * (2.0 * pi * self.sigma_x * self.sigma_y) * fx * fy

    def q_medio(self):
        """Flujo medio sobre la placa, potencia entre área [W/m2]."""
        return self.potencia_total() / self.area

    def factor_interceptacion(self):
        """Fracción de la campana que cae dentro de la placa [-].

        1 - este valor es el derrame (spillage): la parte del cono de radiación
        que se escapa por los bordes y no llega al receptor.
        """
        xc, yc = self.centro
        return (_fraccion_gaussiana(0.0, self.W, xc, self.sigma_x)
                * _fraccion_gaussiana(0.0, self.L, yc, self.sigma_y))

    def factor_pico(self):
        """Relación entre el flujo de pico y el medio, q_pico/q_medio [-].

        Mide lo desigual que es el mapa. Vale 1 si el flujo fuera uniforme y
        crece al estrechar la campana; es el número que gobierna el punto
        caliente del receptor.
        """
        return self.q_pico / self.q_medio()

    # -- Discretización ------------------------------------------------------

    def mapa_nodal(self, nx, ny):
        """Flujo MEDIO de cada celda de una malla nx (anchura) por ny (altura).

        Devuelve una lista de ny filas por nx columnas, ordenada de abajo arriba:
        resultado[j][i] es la celda de la columna i y la fila j, con i creciendo
        en x y j creciendo en y, igual que las coordenadas.

        Cada celda se INTEGRA, no se muestrea en su centro, de modo que
        sum(q_celda * area_celda) reproduce potencia_total() exactamente. Ver la
        nota del módulo.
        """
        if nx < 1 or ny < 1:
            raise ValueError("La malla necesita al menos una celda en cada dirección.")
        xc, yc = self.centro
        dx, dy = self.W / nx, self.L / ny
        area_celda = dx * dy
        # Fracción de campana que cae en cada franja, en cada dirección
        fx = [_fraccion_gaussiana(i * dx, (i + 1) * dx, xc, self.sigma_x) for i in range(nx)]
        fy = [_fraccion_gaussiana(j * dy, (j + 1) * dy, yc, self.sigma_y) for j in range(ny)]
        # La gaussiana 2D es separable: la potencia de la celda es el producto
        # de las dos fracciones por la potencia de la campana entera.
        potencia_campana = self.q_pico * 2.0 * pi * self.sigma_x * self.sigma_y
        return [[potencia_campana * fx[i] * fy[j] / area_celda for i in range(nx)]
                for j in range(ny)]

    def perfil_y(self, ny):
        """Flujo medio de cada una de las ny franjas horizontales [W/m2].

        Es el mapa nodal con una sola columna: el perfil que ve un canal que
        recorre la placa de abajo arriba abarcando toda su anchura.
        """
        return [fila[0] for fila in self.mapa_nodal(1, ny)]

    def resumen(self):
        """Línea de texto con lo esencial del mapa, para los listados."""
        return (
            f"Placa {self.L:.3g} x {self.W:.3g} m | "
            f"sigma = ({self.sigma_x:.3g}, {self.sigma_y:.3g}) m | "
            f"q_pico = {self.q_pico/1e3:.1f} kW/m2 | "
            f"q_medio = {self.q_medio()/1e3:.1f} kW/m2 | "
            f"Q = {self.potencia_total()/1e3:.1f} kW | "
            f"interceptacion = {self.factor_interceptacion()*100:.1f} %"
        )


# =============================================================================
# Mapa por defecto
# =============================================================================

# Placa de 1.5 m de alto por 1.0 m de ancho con 500 kW de pico, un orden de
# magnitud habitual en receptores de torre, y sigmas iguales a un tercio del
# lado correspondiente: la campana queda contenida en la placa con un derrame
# pequeño. Es solo un caso de ejemplo para poder ejecutar el módulo; el mapa de
# verdad lo fijará el dimensionado del campo.
MAPA_EJEMPLO = MapaGaussiano.desde_pico(
    L=1.5, W=1.0, q_pico=500e3, sigma_x=1.0 / 3.0, sigma_y=1.5 / 3.0,
)


# =============================================================================
# Representación y comprobaciones (solo al ejecutar el módulo)
# =============================================================================


def mapa_ascii(mapa, nx=25, ny=15):
    """Dibuja el mapa en la consola, sin dependencias, para verlo de un vistazo."""
    niveles = " .:-=+*#%@"
    celdas = mapa.mapa_nodal(nx, ny)
    q_max = max(max(fila) for fila in celdas)
    print(f"\nMapa de flujo ({nx} x {ny} celdas), de 0 a {q_max/1e3:.0f} kW/m2")
    for fila in reversed(celdas):     # De arriba abajo, como se ve la placa
        linea = "".join(niveles[min(int(q / q_max * len(niveles)), len(niveles) - 1)]
                        for q in fila)
        print("  |" + linea + "|")
    print("  +" + "-" * nx + "+")
    print(f"   x = 0 {' ' * max(nx - 14, 0)} x = W = {mapa.W:.3g} m")


def dibujar(mapa, n=200, archivo=None, mostrar=True, cmap=CMAP, vmin=None, vmax=None):
    """Mapa de color con matplotlib. Devuelve la figura, o None si no lo hay.

    Escala continua, sin líneas de nivel ni bandas: se pinta con imshow sobre
    una malla regular de n x n puntos, no con contourf, para que el degradado
    sea liso y no aparezcan escalones donde no los hay.

    Con archivo se guarda además en disco, que es lo que interesa para meterlo
    en la memoria del TFG. Con mostrar=False no abre ventana, útil para generar
    figuras en lote sin que cada una bloquee la ejecución.

    Por defecto el negro es el mínimo QUE HAY SOBRE LA PLACA, no el flujo nulo.
    Para comparar varias figuras entre sí hay que anclar la escala pasándoles a
    todas el mismo vmin y vmax, en kW/m2: vmin=0 pone el negro en flujo nulo.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no esta instalado: se omite la figura.")
        return None

    xs = [mapa.W * i / (n - 1) for i in range(n)]
    ys = [mapa.L * j / (n - 1) for j in range(n)]
    campo = [[mapa.q(x, y) / 1e3 for x in xs] for y in ys]

    fig, ax = plt.subplots(figsize=(4.0 * mapa.W / mapa.L + 1.5, 4.0))
    imagen = ax.imshow(
        campo, cmap=cmap, origin="lower", aspect="equal",
        extent=(0.0, mapa.W, 0.0, mapa.L), interpolation="bilinear",
        vmin=vmin, vmax=vmax,
    )
    ax.set_xlabel("x, anchura [m]")
    ax.set_ylabel("y, altura [m]")
    ax.set_title("Flujo incidente [kW/m$^2$]")
    fig.colorbar(imagen, ax=ax)
    fig.tight_layout()
    if archivo:
        fig.savefig(archivo, dpi=200)
        print(f"Figura guardada en {archivo}")
    if mostrar:
        plt.show()      # Bloquea hasta que se cierra la ventana
    return fig


def comprobar(mapa=MAPA_EJEMPLO):
    """Comprueba que el mapa conserva la energía y que las integrales cuadran."""
    Q = mapa.potencia_total()
    print(f"\n{mapa.resumen()}")
    xc, yc = mapa.centro
    print(f"q en el centro {mapa.q(xc, yc)/1e3:10.2f} kW/m2 (q_pico "
          f"{mapa.q_pico/1e3:.2f})")
    print(f"q en una esquina{mapa.q(0.0, 0.0)/1e3:10.2f} kW/m2")

    encabezado = (f"\n{'malla':>12}{'suma nodal [kW]':>18}{'error [%]':>12}"
                  f"{'q_max nodal':>14}{'muestreo [%]':>14}")
    print(encabezado)
    print('-' * len(encabezado.strip('\n')))
    for nx, ny in ((1, 1), (2, 3), (5, 8), (10, 15), (40, 60), (200, 300)):
        celdas = mapa.mapa_nodal(nx, ny)
        area_celda = mapa.area / (nx * ny)
        suma = sum(q for fila in celdas for q in fila) * area_celda
        # Lo que daría muestrear el centro de cada celda en vez de integrarla
        muestreo = sum(
            mapa.q(mapa.W * (i + 0.5) / nx, mapa.L * (j + 0.5) / ny)
            for j in range(ny) for i in range(nx)
        ) * area_celda
        q_max = max(max(fila) for fila in celdas)
        print(f"{f'{nx} x {ny}':>12}{suma/1e3:18.4f}{(suma/Q - 1)*100:12.2e}"
              f"{q_max/1e3:14.2f}{(muestreo/Q - 1)*100:+14.2f}")

    # Contraste con una cuadratura numérica independiente, si hay scipy
    try:
        from scipy.integrate import dblquad
    except ImportError:
        print("\nscipy no esta instalado: se omite la cuadratura de contraste.")
    else:
        numerica, _ = dblquad(lambda y, x: mapa.q(x, y), 0.0, mapa.W, 0.0, mapa.L)
        print(f"\nPotencia por cuadratura numerica: {numerica/1e3:.4f} kW "
              f"(erf: {Q/1e3:.4f} kW, error {(Q/numerica - 1)*100:+.2e} %)")

    # El error clásico: integrar la campana en todo el plano en vez de en la placa
    infinita = mapa.q_pico * 2.0 * pi * mapa.sigma_x * mapa.sigma_y
    print(f"Campana completa (plano infinito): {infinita/1e3:.4f} kW, "
          f"un {(infinita/Q - 1)*100:+.1f} % sobre la placa")

    # Ida y vuelta: fijar la potencia y recuperar el mismo pico
    rehecho = MapaGaussiano.desde_potencia(mapa.L, mapa.W, Q, mapa.sigma_x, mapa.sigma_y)
    print(f"desde_potencia({Q/1e3:.2f} kW) devuelve q_pico = "
          f"{rehecho.q_pico/1e3:.4f} kW/m2 (original {mapa.q_pico/1e3:.4f})")


if __name__ == "__main__":
    comprobar()

    # Perfil que ve un canal que recorre la placa de abajo arriba
    print("\nPerfil en altura, 10 franjas [kW/m2]:")
    print("  " + "  ".join(f"{q/1e3:6.1f}" for q in MAPA_EJEMPLO.perfil_y(10)))

    # Campana ancha: aquí el derrame deja de ser despreciable
    ancho = MapaGaussiano.desde_pico(L=1.5, W=1.0, q_pico=500e3, sigma_x=0.8, sigma_y=1.2)
    print(f"\nCampana ancha: {ancho.resumen()}")
    infinita = ancho.q_pico * 2.0 * pi * ancho.sigma_x * ancho.sigma_y
    print(f"  integrar en el plano infinito daria {infinita/1e3:.1f} kW, "
          f"un {(infinita/ancho.potencia_total() - 1)*100:+.0f} % de mas")

    # La figura, al final: plt.show() bloquea hasta que se cierra la ventana, y
    # así todo lo anterior ya está impreso. Si no hay matplotlib se recurre al
    # mapa de consola, que para eso se escribió (y sirve también por SSH, sin
    # entorno gráfico).
    if dibujar(MAPA_EJEMPLO) is None:
        mapa_ascii(MAPA_EJEMPLO)
