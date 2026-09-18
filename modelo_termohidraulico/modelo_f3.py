"""Canal interno de un PPHE, con la velocidad obtenida del caudal (ec. 3).

Modelo directo, sin iterar: se fijan el número de chapas y la geometría del
panel, y la velocidad sale de la ecuación (3) de O. Arsenyeva et al.,

    w = G / (rho * N * f_chI)

de donde salen Re, Pr, Nu, el factor de fricción y el coeficiente de película.
La pérdida de carga es un RESULTADO, obtenido de la ec. (19); no se impone.

Es la alternativa a modelo_iterativo.py, que hace lo contrario: impone la pérdida
de carga admisible y despeja la velocidad de la ecuación implícita (28). Las dos
vías son equivalentes y dan el mismo punto de diseño cuando se cierran con las
mismas condiciones; ésta tiene la ventaja de partir del caudal, que es el dato
duro en un receptor solar, y de no depender de la constante discutida de la (28).

CADA CASO LLEVA SU PROPIO FLUIDO. Es deliberado: la geometría y las condiciones
de proceso proceden de casos de estudio distintos del artículo, y tenerlas
separadas de la lista de casos hacía fácil emparejar la geometría de uno con el
fluido de otro sin que nada avisara.

NOTAS SOBRE LOS DATOS DE ENTRADA:

  - La columna n_pl de las tablas del artículo es el número de CHAPAS, no de
    canales. Cada panel se forma soldando dos chapas, y por eso la ec. (13)
    obliga a que n_pl sea par; los canales internos son n_pl/2. Ver
    CHAPAS_POR_PANEL para probar la otra lectura.

  - La separación entre paneles b NO interviene en ningún cálculo de este
    modelo: tanto deI (ec. 14) como f_chI (ec. 16) dependen solo de b_i. Se
    conserva en cada caso porque define la geometría del canal externo, que aquí
    no se calcula, y porque es la variable de diseño del artículo.

  - La longitud L_F se toma de la tabla de referencia. Sin canal externo no hay
    coeficiente global U y, por tanto, la ec. (26) no puede darla. Solo se usa
    para la pérdida de carga: no afecta a la velocidad ni a ningún adimensional.
"""

from dataclasses import dataclass, replace
import formulas as form
import asignacion as asig
import geometria as geom

CHAPAS_POR_PANEL = 2      # Canales internos = n_pl / CHAPAS_POR_PANEL; a 1 para la lectura alternativa


# =============================================================================
# Fluidos
# =============================================================================

# Caso de estudio 1 (Tabla 4). El canal interno lleva el fluido frío.
AGUA_FRIA = dict(
    nombre="Agua (frío)", rho=980, mu=0.8e-3, k=0.618, cp=4175,
    t_in=10, t_out=50, G=30,
)

# Caso de estudio 2 (Tabla 7). Aquí el canal interno lleva el CALIENTE, que es
# el de menor caudal. Propiedades tomadas de la ref. [24] del artículo.
PRODUCTO_1 = dict(
    nombre="Producto 1 (caliente)", rho=780, mu=0.04443e-3, k=0.0714, cp=2900,
    t_in=177.4, t_out=155.6, G=37.0/3.6,
)


# =============================================================================
# Casos
# =============================================================================

# n_pl, b y L_F son los del diseño óptimo que publica el artículo, y el bloque
# ref recoge lo que reporta para ese mismo canal, para contrastar. La Tabla 9 no
# publica Re ni h por corriente, de ahí los None.
CASOS = [
    dict(nombre="PPHE1", geometria=geom.PPHE1, fluido=AGUA_FRIA, tabla="Tabla 5",
         n_pl=70, b=5.5e-3, L_F=2.072,
         ref=dict(w=1.835, Re=1.081e4, del_P=6.000e4, h=9804)),
    dict(nombre="PPHE2", geometria=geom.PPHE2, fluido=AGUA_FRIA, tabla="Tabla 5",
         n_pl=112, b=7.5e-3, L_F=1.200,
         ref=dict(w=1.280, Re=6.670e3, del_P=6.000e4, h=11220)),
    dict(nombre="PPHE3", geometria=geom.PPHE3, fluido=AGUA_FRIA, tabla="Tabla 5",
         n_pl=66, b=20e-3, L_F=2.230,
         ref=dict(w=0.940, Re=1.139e4, del_P=6.000e4, h=10206)),
    dict(nombre="HE1", geometria=geom.PPHE2, fluido=PRODUCTO_1, tabla="Tabla 9",
         n_pl=24, b=1.5e-3, L_F=2.250,
         ref=dict(w=1.943, Re=None, del_P=1.000e5, h=None)),
]


# =============================================================================
# Estructuras de datos
# =============================================================================


@dataclass(frozen=True)
class PPHEResult:
    nombre: str
    tabla: str
    fluido: str
    N: float
    f_ch_int: float
    d_h_int: float
    u_m_int: float
    Pr_int: float
    Re_int: float
    Nu_int: float
    f_int: float
    h_int: float
    del_P_int: float
    N_implicito: float
    ref: dict


def _fmt(valor, ancho, decimales):
    """Formatea un valor de referencia, o un guion si el artículo no lo publica."""
    return f"{'-':>{ancho}}" if valor is None else f"{valor:>{ancho}.{decimales}f}"


def imprimir_resultados(resultados):
    """Da formato a la impresión por consola de los resultados."""
    encabezado = (
        f"{'Caso':<8}{'N':>6}{'f_ch [m2]':>13}{'Dh [m]':>14}{'u [m/s]':>11}"
        f"{'Pr':>9}{'Re':>12}{'Nu':>11}{'f':>10}{'h [W/m2K]':>13}{'dP [kPa]':>12}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for r in resultados:
        print(
            f"{r.nombre:<8}{r.N:6.0f}{r.f_ch_int:13.6e}{r.d_h_int:14.6e}{r.u_m_int:11.4f}"
            f"{r.Pr_int:9.4f}{r.Re_int:12.4g}{r.Nu_int:11.4f}{r.f_int:10.4f}"
            f"{r.h_int:13.2f}{r.del_P_int/1e3:12.2f}"
        )


def imprimir_comparacion(resultados):
    """Contrasta lo calculado con lo que reporta la tabla de referencia."""
    encabezado = (
        f"{'Caso':<8}{'Referencia':<12}{'u [m/s]':>10}{'ref':>9}{'ratio':>8}"
        f"{'Re':>11}{'ref':>11}{'h [W/m2K]':>12}{'ref':>10}"
        f"{'dP [kPa]':>11}{'ref':>9}{'N':>7}{'N implic.':>11}{'ratio':>8}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for r in resultados:
        ref = r.ref
        print(
            f"{r.nombre:<8}{r.tabla:<12}{r.u_m_int:10.4f}{_fmt(ref['w'], 9, 3)}"
            f"{ref['w']/r.u_m_int:8.3f}"
            f"{r.Re_int:11.4g}{_fmt(ref['Re'], 11, 0)}"
            f"{r.h_int:12.2f}{_fmt(ref['h'], 10, 0)}"
            f"{r.del_P_int/1e3:11.2f}{_fmt(ref['del_P']/1e3, 9, 1)}"
            f"{r.N:7.0f}{r.N_implicito:11.2f}{r.N/r.N_implicito:8.3f}"
        )


# =============================================================================
# Cálculos numéricos
# =============================================================================

resultados = []

for caso in CASOS:
    nombre = caso["nombre"]
    fluido = caso["fluido"]

    # b del diseño de referencia, por si el fichero de geometría llevara otro
    i = replace(caso["geometria"], b=caso["b"])

    constantes_n = asig.asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i)

    """VELOCIDAD (ecuación (3): el caudal repartido entre N canales)"""

    N = caso["n_pl"]/CHAPAS_POR_PANEL                                                                # Canales, no chapas
    f_ch_int = form.f_chI(i.b_i, i.w_pp, i.w_e)                                                      # Sección de paso (16)
    u_m_int = form.velocidad_ec3(G=fluido["G"], rho=fluido["rho"], N=N, f_ch=f_ch_int)               # Velocidad media (3)

    """CÁLCULO INTERNO"""

    d_h_int = form.deI(i.b_i)                                                                        # Diámetro equivalente (14)
    Re_int = form.Re(rho=fluido["rho"], u=u_m_int, dh=d_h_int, mu=fluido["mu"])                       # Número de Reynolds
    Pr_int = form.Pr(cp=fluido["cp"], mu=fluido["mu"], k=fluido["k"])                                 # Número de Prandtl
    Nu_int = form.NuI(n3=constantes_n.n3, n4=constantes_n.n4, n5=constantes_n.n5,
                      Re=Re_int, Pr=Pr_int)                                                          # Número de Nusselt (20)
    f_int = form.fI(n1=constantes_n.n1, Re=Re_int, n2=constantes_n.n2)                                # Factor de fricción (18)
    h_int = Nu_int * fluido["k"] / d_h_int                                                            # Coef. de transferencia de calor

    """PÉRDIDA DE CARGA (resultado, no dato: L_F es la de la tabla de referencia)"""

    del_P_int = form.perdida_de_carga(f_int, caso["L_F"], d_h_int, fluido["rho"], u_m_int)            # Pérdida de carga (19)

    # Canales que harían falta para dar la velocidad que reporta el artículo
    N_implicito = fluido["G"]/(fluido["rho"]*caso["ref"]["w"]*f_ch_int)

    print(f"\n--- {nombre} ({caso['tabla']}) ---")
    print(f"{fluido['nombre']}, G={fluido['G']:.4f} kg/s")
    print(f"n_pl={caso['n_pl']} chapas -> N={N:.0f} canales, b={caso['b']*1e3:.1f} mm (no interviene en el canal interno)")
    print(f"u={u_m_int:.4f} m/s (referencia: {caso['ref']['w']:.3f}, desviación {u_m_int/caso['ref']['w']*100-100:+.2f} %)")
    print(f"dP={del_P_int/1e3:.2f} kPa con L_F={caso['L_F']} m (referencia: {caso['ref']['del_P']/1e3:.0f} kPa)")

    resultados.append(
        PPHEResult(
            nombre=nombre,
            tabla=caso["tabla"],
            fluido=fluido["nombre"],
            N=N,
            f_ch_int=f_ch_int,
            d_h_int=d_h_int,
            u_m_int=u_m_int,
            Pr_int=Pr_int,
            Re_int=Re_int,
            Nu_int=Nu_int,
            f_int=f_int,
            h_int=h_int,
            del_P_int=del_P_int,
            N_implicito=N_implicito,
            ref=caso["ref"],
        )
    )

imprimir_resultados(resultados)
imprimir_comparacion(resultados)
