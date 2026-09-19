"""Propiedades termofísicas del fluido de trabajo del receptor, en función de T.

Cada fluido es un objeto que sabe devolver su densidad, calor específico,
conductividad térmica y viscosidad dinámica a la temperatura (y presión) que se
le pida. De momento solo está el AIRE; para añadir otro gas basta con crear otro
GasIdealPolinomico con sus coeficientes, sin tocar el resto del modelo.

UNIDADES. Todo en SI y las TEMPERATURAS EN KELVIN, sin excepción. El resto de
módulos del TFG manejan grados Celsius en las condiciones de proceso, así que al
llamar aquí hay que sumar CERO_CELSIUS. Se ha preferido K a C porque la ley de
gases ideales y los ajustes polinómicos son en K, y convertir dentro de cada
correlación invitaba a olvidarse de hacerlo en alguna.

    rho  densidad                  [kg/m3]
    cp   calor específico a p cte  [J/(kg*K)]
    k    conductividad térmica     [W/(m*K)]
    mu   viscosidad dinámica       [Pa*s]

DE DÓNDE SALEN LOS NÚMEROS. La densidad es gas ideal, rho = p/(R_esp*T). Las
otras tres son ajustes polinómicos de grado 5 en x = T/1000 a los datos de
referencia de CoolProp 7.2.0 para aire seco, que a su vez implementa

    E.W. Lemmon, R.T. Jacobsen, S.G. Penoncello, D.G. Friend, "Thermodynamic
    properties of air and mixtures of nitrogen, argon and oxygen from 60 to
    2000 K at pressures to 2000 MPa", J. Phys. Chem. Ref. Data 29 (2000)
    331-385. [ecuación de estado, de donde sale cp]

    E.W. Lemmon, R.T. Jacobsen, "Viscosity and thermal conductivity equations
    for nitrogen, oxygen, argon and air", Int. J. Thermophys. 25 (2004) 21-69.
    [mu y k]

Los ajustes se rehacen y se comprueban ejecutando este módulo (ver el bloque
final): imprime la desviación de cada propiedad frente a CoolProp. En 250-1600 K
es 0.13 % en cp, 0.07 % en mu y 0.05 % en k, dos órdenes de magnitud por debajo
de la incertidumbre de las correlaciones de Nu y f del modelo termohidráulico.

POR QUÉ UN AJUSTE Y NO LLAMAR A CoolProp EN CADA NODO. Porque el receptor se
resuelve nodo a nodo e iterando, y cada evaluación de CoolProp cuesta del orden
de cien veces más que un polinomio de Horner. CoolProp queda como referencia de
validación, no como dependencia de ejecución: el modelo corre sin tenerlo
instalado.

POR QUÉ NO SUTHERLAND. La ley de Sutherland con las constantes clásicas
(mu0 = 1.716e-5 Pa*s, S = 110.4 K) es cómoda, pero se desvía hasta un 6.6 % en
mu y un 6.8 % en k a 1500 K, y el receptor trabaja justo en esa zona. Reajustar
sus dos constantes al rango no lo arregla (quedan 3.6 % y 7.8 %): la forma
funcional no da más de sí a alta temperatura. De ahí el polinomio.

VALIDEZ EN PRESIÓN. cp, k y mu son los del límite de gas diluido, es decir, los
de 1 atm: se desprecia su dependencia con la presión. Comparando con CoolProp,
el error que eso introduce es <=1.4 % hasta 10 bar (el caso peor es cp a 300 K y
10 bar) y llega al 3 % a 20 bar. La densidad de gas ideal se desvía <0.5 % hasta
20 bar. Si el receptor se presuriza por encima de unos 10 bar habrá que
sustituir estas correlaciones por CoolProp o añadirles el término de densidad.
"""

from dataclasses import dataclass
from typing import Sequence
import warnings

# =============================================================================
# Constantes
# =============================================================================

R_UNIVERSAL = 8.314462618     # Constante universal de los gases [J/(mol*K)] (CODATA 2018)
CERO_CELSIUS = 273.15         # Conversión C -> K [K]
P_ATM = 101325.0              # Presión atmosférica normal, presión por defecto [Pa]
T_REF_ENTALPIA = 273.15       # Temperatura de referencia de la entalpía, h(T_REF)=0 [K]

# Qué hacer si se pide una temperatura fuera del rango de ajuste del fluido:
# False -> ValueError (por defecto; un polinomio de grado 5 extrapolado se
# dispara muy deprisa). True -> solo avisa y extrapola, útil para depurar una
# iteración que se sale de rango en los primeros pasos.
EXTRAPOLAR = False


# =============================================================================
# Estado termofísico
# =============================================================================


@dataclass(frozen=True)
class Propiedades:
    """Las cuatro propiedades y sus derivadas, en un punto (T, p) concreto.

    Se devuelven juntas para no repetir el chequeo de rango ni la evaluación de
    los polinomios en cada nodo.
    """

    fluido: str
    T: float        # Temperatura [K]
    p: float        # Presión [Pa]
    rho: float      # Densidad [kg/m3]
    cp: float       # Calor específico a presión constante [J/(kg*K)]
    k: float        # Conductividad térmica [W/(m*K)]
    mu: float       # Viscosidad dinámica [Pa*s]

    @property
    def nu(self):
        """Viscosidad cinemática, mu/rho [m2/s]."""
        return self.mu / self.rho

    @property
    def alfa(self):
        """Difusividad térmica, k/(rho*cp) [m2/s]."""
        return self.k / (self.rho * self.cp)

    @property
    def Pr(self):
        """Número de Prandtl, cp*mu/k [-]."""
        return self.cp * self.mu / self.k

    @property
    def t_celsius(self):
        """La misma temperatura, en grados Celsius [C]."""
        return self.T - CERO_CELSIUS


# =============================================================================
# Fluido genérico: gas ideal con propiedades de transporte polinómicas
# =============================================================================


def _horner(coefs: Sequence[float], x: float) -> float:
    """Evalúa un polinomio dado por sus coeficientes ASCENDENTES (c0 + c1*x + ...)."""
    resultado = 0.0
    for c in reversed(coefs):
        resultado = resultado * x + c
    return resultado


@dataclass(frozen=True)
class GasIdealPolinomico:
    """Gas ideal cuyas cp, k y mu son polinomios en x = T/1000.

    Los coeficientes van en orden ASCENDENTE y en las unidades que indica cada
    campo; el escalado a SI lo hace el método correspondiente. Se trabaja con
    x = T/1000 y no con T porque así los coeficientes quedan del orden de la
    unidad y el ajuste no pierde cifras significativas.
    """

    nombre: str
    M: float                    # Masa molar [kg/mol]
    coef_cp: Sequence[float]    # cp [J/(kg*K)] = poly(x)
    coef_k: Sequence[float]     # k  [mW/(m*K)] = poly(x)
    coef_mu: Sequence[float]    # mu [uPa*s]    = poly(x)
    T_min: float                # Límite inferior del ajuste [K]
    T_max: float                # Límite superior del ajuste [K]
    fuente: str = ""            # De dónde salen los datos ajustados

    # -- Comprobaciones ------------------------------------------------------

    def _verificar(self, T):
        """Impide extrapolar el ajuste fuera del rango en que se hizo.

        Ver EXTRAPOLAR para relajarlo.
        """
        if self.T_min <= T <= self.T_max:
            return
        mensaje = (
            f"{self.nombre}: T = {T:.2f} K = {T - CERO_CELSIUS:.2f} C queda fuera "
            f"del rango de ajuste [{self.T_min:.0f}, {self.T_max:.0f}] K."
        )
        if EXTRAPOLAR:
            warnings.warn(mensaje + " Se extrapola.", stacklevel=3)
        else:
            raise ValueError(mensaje)

    @property
    def R(self):
        """Constante específica del gas, R_universal/M [J/(kg*K)]."""
        return R_UNIVERSAL / self.M

    # -- Las cuatro propiedades ----------------------------------------------

    def rho(self, T, p=P_ATM):
        """Densidad por la ley de los gases ideales, rho = p/(R*T) [kg/m3]."""
        self._verificar(T)
        return p / (self.R * T)

    def cp(self, T):
        """Calor específico a presión constante [J/(kg*K)]. Ajuste, ver el módulo."""
        self._verificar(T)
        return _horner(self.coef_cp, T / 1000.0)

    def k(self, T):
        """Conductividad térmica [W/(m*K)]. Ajuste, ver el módulo."""
        self._verificar(T)
        return _horner(self.coef_k, T / 1000.0) * 1e-3

    def mu(self, T):
        """Viscosidad dinámica [Pa*s]. Ajuste, ver el módulo."""
        self._verificar(T)
        return _horner(self.coef_mu, T / 1000.0) * 1e-6

    # -- Derivadas -----------------------------------------------------------

    def Pr(self, T):
        """Número de Prandtl, cp*mu/k [-]."""
        return self.cp(T) * self.mu(T) / self.k(T)

    def nu(self, T, p=P_ATM):
        """Viscosidad cinemática, mu/rho [m2/s]."""
        return self.mu(T) / self.rho(T, p)

    def alfa(self, T, p=P_ATM):
        """Difusividad térmica, k/(rho*cp) [m2/s]."""
        return self.k(T) / (self.rho(T, p) * self.cp(T))

    def propiedades(self, T, p=P_ATM):
        """Las cuatro propiedades de una vez, en un objeto Propiedades."""
        self._verificar(T)
        x = T / 1000.0
        return Propiedades(
            fluido=self.nombre,
            T=T,
            p=p,
            rho=p / (self.R * T),
            cp=_horner(self.coef_cp, x),
            k=_horner(self.coef_k, x) * 1e-3,
            mu=_horner(self.coef_mu, x) * 1e-6,
        )

    # -- Entalpía ------------------------------------------------------------

    def h(self, T):
        """Entalpía específica referida a T_REF_ENTALPIA [J/kg].

        Es la integral analítica del polinomio de cp, no cp*(T-T_ref): entre 300
        y 1200 K el cp del aire sube un 17 %, así que cerrar el balance con el cp
        de la entrada se deja un 7.4 % del calor absorbido por el camino.

        Como el origen es arbitrario, lo que tiene sentido físico es siempre una
        DIFERENCIA de entalpías, h(T2)-h(T1).
        """
        self._verificar(T)
        return self._integral_cp(T / 1000.0) - self._integral_cp(T_REF_ENTALPIA / 1000.0)

    def _integral_cp(self, x):
        """Primitiva del polinomio de cp respecto de T, evaluada en x = T/1000.

        El factor 1000 es el cambio de variable dT = 1000*dx.
        """
        return 1000.0 * sum(c * x ** (i + 1) / (i + 1) for i, c in enumerate(self.coef_cp))

    def cp_medio(self, T1, T2):
        """Calor específico medio entre T1 y T2 [J/(kg*K)].

        Es el cp que hay que usar en Q = G*cp_medio*(T2-T1) para que el calor
        salga exacto, y el que cierra el balance de un tramo con salto térmico
        grande.
        """
        if abs(T2 - T1) < 1e-9:
            return self.cp(T1)
        return (self.h(T2) - self.h(T1)) / (T2 - T1)


# =============================================================================
# Fluidos disponibles
# =============================================================================

# Aire seco. Coeficientes ajustados a CoolProp 7.2.0 ('Air') a 1 atm entre 250 y
# 1600 K, por mínimos cuadrados con peso 1/y, es decir, minimizando el error
# RELATIVO y no el absoluto: las tres propiedades varían por un factor ~3 en el
# rango, y un ajuste sin pesar descuidaría el extremo frío.
#
# Desviación máxima frente a CoolProp en 250-1600 K: cp 0.13 %, mu 0.07 %, k 0.05 %.
AIRE = GasIdealPolinomico(
    nombre="Aire",
    M=28.96546e-3,
    coef_cp=(1086.4453, -659.81616, 1702.3009, -1486.6131, 586.96592, -88.252623),
    coef_k=(-0.25684666, 107.47369, -77.64382, 58.104571, -24.214681, 4.2141488),
    coef_mu=(0.45586184, 75.791245, -64.487066, 48.07116, -20.037935, 3.4861614),
    T_min=250.0,
    T_max=1600.0,
    fuente="Lemmon et al. (2000, 2004), vía CoolProp 7.2.0",
)

FLUIDOS = {"aire": AIRE}


def fluido(nombre="aire"):
    """Devuelve el fluido de trabajo a partir de su nombre."""
    try:
        return FLUIDOS[nombre.lower()]
    except KeyError:
        raise ValueError(
            f"Fluido '{nombre}' no disponible. Hay: {', '.join(sorted(FLUIDOS))}."
        ) from None


# =============================================================================
# Tabla y validación (solo al ejecutar el módulo)
# =============================================================================


def imprimir_tabla(f=AIRE, T_ini=300.0, T_fin=1500.0, paso=100.0, p=P_ATM):
    """Tabla de propiedades del fluido, para mirarla de un vistazo."""
    print(f"\n{f.nombre} a {p/1e5:.3f} bar   (R = {f.R:.3f} J/(kg*K))")
    encabezado = (
        f"{'T [K]':>8}{'t [C]':>8}{'rho':>10}{'cp':>9}{'k':>9}{'mu':>12}"
        f"{'nu':>12}{'alfa':>12}{'Pr':>8}"
    )
    print(encabezado)
    print('-' * len(encabezado))
    T = T_ini
    while T <= T_fin + 1e-9:
        pr = f.propiedades(T, p)
        print(
            f"{pr.T:8.0f}{pr.t_celsius:8.1f}{pr.rho:10.4f}{pr.cp:9.1f}{pr.k:9.5f}"
            f"{pr.mu:12.4e}{pr.nu:12.4e}{pr.alfa:12.4e}{pr.Pr:8.4f}"
        )
        T += paso


def validar_con_coolprop(f=AIRE, p=P_ATM, n=271):
    """Contrasta el ajuste con CoolProp, si está instalado.

    Recorre todo el rango de validez del fluido y da la desviación máxima de
    cada propiedad. Es la comprobación que respalda los porcentajes citados en
    el docstring del módulo: si algún día se cambian los coeficientes, basta con
    volver a ejecutar este módulo para saber qué error tienen.
    """
    try:
        from CoolProp.CoolProp import PropsSI
    except ImportError:
        print("\nCoolProp no esta instalado: se omite la validacion.")
        return None

    temperaturas = [f.T_min + (f.T_max - f.T_min) * i / (n - 1) for i in range(n)]
    errores = {"rho": [], "cp": [], "k": [], "mu": []}
    for T in temperaturas:
        pr = f.propiedades(T, p)
        ref = dict(
            rho=PropsSI('D', 'T', T, 'P', p, 'Air'),
            cp=PropsSI('C', 'T', T, 'P', p, 'Air'),
            k=PropsSI('L', 'T', T, 'P', p, 'Air'),
            mu=PropsSI('V', 'T', T, 'P', p, 'Air'),
        )
        for nombre, valor in ref.items():
            errores[nombre].append((getattr(pr, nombre) / valor - 1) * 100)

    print(f"\nValidacion frente a CoolProp en [{f.T_min:.0f}, {f.T_max:.0f}] K "
          f"a {p/1e5:.3f} bar, {n} puntos")
    encabezado = f"{'Propiedad':<12}{'err max [%]':>13}{'en T [K]':>11}{'err medio [%]':>15}"
    print(encabezado)
    print('-' * len(encabezado))
    for nombre, e in errores.items():
        i = max(range(len(e)), key=lambda j: abs(e[j]))
        medio = sum(abs(v) for v in e) / len(e)
        print(f"{nombre:<12}{e[i]:13.4f}{temperaturas[i]:11.0f}{medio:15.4f}")
    return errores


if __name__ == "__main__":
    imprimir_tabla()
    validar_con_coolprop()

    # Por qué hace falta cp_medio: con el salto térmico de un receptor, tomar el
    # cp de la entrada se queda corto en el calor absorbido en casi un 8 %. El cp
    # de la temperatura media, en cambio, sí es una buena aproximación del medio
    # (el error es el término (1/24)*cp''*dT**2 del desarrollo, despreciable
    # aquí); cp_medio lo hace exacto y de paso vale para cualquier salto.
    T1, T2 = 300.0, 1200.0
    cp_m = AIRE.cp_medio(T1, T2)
    print(f"\nSalto de {T1:.0f} a {T2:.0f} K:")
    print(f"  cp medio (integrado)        {cp_m:8.1f} J/(kg*K)")
    print(f"  cp a la T media ({(T1+T2)/2:.0f} K)      {AIRE.cp((T1+T2)/2):8.1f} J/(kg*K)")
    print(f"  cp a la entrada ({T1:.0f} K)       {AIRE.cp(T1):8.1f} J/(kg*K)"
          f"  -> {AIRE.cp(T1)/cp_m*100-100:+.1f} % en el calor")
    print(f"  salto de entalpia           {(AIRE.h(T2) - AIRE.h(T1))/1e3:8.1f} kJ/kg")
