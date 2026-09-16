"""Canal interno de un PPHE con la velocidad obtenida del caudal (ec. 3).

Réplica del PRIMER caso de estudio del artículo de O. Arsenyeva et al. (agua,
Tabla 4), contrastada contra la Tabla 5. Se calcula únicamente el canal INTERNO,
por el que circula el fluido frío.

No hay método iterativo y no se impone ninguna pérdida de carga. Los datos de
entrada son el número de chapas y la separación entre paneles de la Tabla 5, y
la velocidad sale directamente de la ecuación (3):

    w = G / (rho * N * f_chI)

de donde salen Re, Pr, Nu, el factor de fricción y el coeficiente de película.
La pérdida de carga es un RESULTADO, obtenido de la ec. (19).

Notas sobre los datos de entrada:

  - La columna n_pl de la Tabla 5 es el número de CHAPAS, no de canales. Cada
    panel se forma soldando dos chapas, y de ahí que la ec. (13) obligue a que
    n_pl sea par, así que los canales internos son n_pl/2. La constante
    CHAPAS_POR_PANEL permite probar la otra lectura.

  - La separación entre paneles b no interviene en ningún cálculo del canal
    interno: tanto d_eI (ec. 14) como f_chI (ec. 16) dependen solo de b_i. Se
    recoge en la tabla de casos por completitud, y porque fija la geometría del
    canal externo, que aquí no se calcula.

  - La longitud L_F se toma de la Tabla 5 porque sin el canal externo no hay
    coeficiente global U y, por tanto, no puede obtenerse de la ec. (26). Solo
    se usa para la pérdida de carga.

Correlaciones tomadas de O. Arsenyeva et al.
"""

from dataclasses import dataclass, replace
import formulas as form
import asignacion as asig
import geometría as geom

CHAPAS_POR_PANEL = 2      # canales internos = n_pl / CHAPAS_POR_PANEL ; ponlo a 1 para la otra lectura

# =============================================================================
# Caso de estudio 1: fluido frío (Tabla 4) y diseños óptimos (Tabla 5)
# =============================================================================

              #   kg/m3       Pa*s      w/(m*K)   J/(kg*K)      °C         °C        kg/s
#props_fl_int = {"rho":980, "mu":0.0008, "k":0.618, "cp":4175, "t_in":10, "t_out":50, "G":30}  #Fluido 2 del artículo: interno, frío
props_fl_int = {"rho":780, "mu":0.04443e-3, "k":0.0714, "cp":2900, "t_in":177.4, "t_out":155.6, "G":10.278} #Producto 2 (Tabla 7)


# n_pl, b y L_F son de la Tabla 5; el bloque t5 recoge lo que el artículo
# reporta del lado frío, para contrastar.
CASOS = [
    dict(nombre="PPHE1", geometria=geom.PPHE1, n_pl=70,  b=5.5e-3, L_F=2.072,
         t5=dict(w=1.835, Re=1.081e4, del_P=6.000e4, h=9804)),
    dict(nombre="PPHE2", geometria=geom.PPHE2, n_pl=24, b=1.5e-3, L_F=2.25,
         t5=dict(w=1.280, Re=6.670e3, del_P=6.000e4, h=11220)),
    dict(nombre="PPHE3", geometria=geom.PPHE3, n_pl=66,  b=20e-3,  L_F=2.230,
         t5=dict(w=0.940, Re=1.139e4, del_P=6.000e4, h=10206)),
]

# =============================================================================
# Estructuras de datos
# =============================================================================



@dataclass(frozen=True)
class PPHEResult:
    nombre: str
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
    t5: dict


def imprimir_resultados(resultados):
    """Da formato a la impresión por consola de los resultados."""
    encabezado = (
        f"{'Caso':<8}{'N':>7}{'f_ch [m2]':>13}{'Dh (in) [m]':>15}{'u (in) [m/s]':>15}{'Pr (in)':>10}{'Re (in)':>15}{'Nu (in)':>12}{'f (in)':>12}{'h (in) [W/m2K]':>16}{'dP (in) [kPa]':>15}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for res in resultados:
        print(
            f"{res.nombre:<8}{res.N:7.0f}{res.f_ch_int:13.6e}{res.d_h_int:15.6e}{res.u_m_int:15.4f}{res.Pr_int:10.4f}{res.Re_int:15.2e}{res.Nu_int:12.4f}{res.f_int:12.4f}{res.h_int:16.2f}{res.del_P_int/1e3:15.2f}"
        )




# =============================================================================
# Cálculos numéricos
# =============================================================================

resultados = []

for caso in CASOS:
    nombre = caso["nombre"]

    # b de la Tabla 5, por si el fichero de geometría llevara otro valor
    i = replace(caso["geometria"], b=caso["b"])

    constantes_n = asig.asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i)

    """VELOCIDAD (ecuación (3): el caudal repartido entre N canales)"""

    N = caso["n_pl"]/CHAPAS_POR_PANEL                                                                                           # Canales, no chapas
    f_ch_int = form.f_chI(i.b_i, i.w_pp, i.w_e)                                                                                 # Sección de paso (16)
    u_m_int = form.velocidad_ec3(G=props_fl_int.get("G"), rho=props_fl_int.get("rho"), N=N, f_ch=f_ch_int)                       # Velocidad media (3)

    """CÁLCULO INTERNO"""

    d_h_int = form.deI(i.b_i)                                                                                                   # Diámetro equivalente
    Re_int = form.Re(rho=props_fl_int.get("rho"), u=u_m_int, dh=d_h_int, mu=props_fl_int.get("mu"))                              # Número de Reynolds
    Pr_int = form.Pr(cp=props_fl_int.get("cp"), mu=props_fl_int.get("mu"), k=props_fl_int.get("k"))                              # Número de Prandtl
    Nu_int = form.NuI(n3=constantes_n[2], n4=constantes_n[3], n5=constantes_n[4], Re=Re_int, Pr=Pr_int)                          # Número de Nusselt
    f_int = form.fI(n1=constantes_n[0], Re=Re_int, n2=constantes_n[1])                                                           # Factor de fricción
    h_int = Nu_int * props_fl_int.get("k") / d_h_int                                                                             # Coef. de transferencia de calor

    """PÉRDIDA DE CARGA (resultado, no dato: L_F es la que reporta la Tabla 5)"""

    del_P_int = form.perdida_de_carga(f_int, caso["L_F"], d_h_int, props_fl_int.get("rho"), u_m_int)

    resultados.append(
        PPHEResult(
            nombre=nombre,
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
            t5=caso["t5"],
        )
    )

imprimir_resultados(resultados)