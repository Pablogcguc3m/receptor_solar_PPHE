"""Modelo termohidráulico de PPHE con la velocidad obtenida del caudal (ec. 3).

Alternativa a modelo_termohidraulico.py, contrastada contra el SEGUNDO caso de
estudio del artículo (tren de precalentamiento de crudo): condiciones de proceso
de la Tabla 7 y geometrías resultantes de la Tabla 9.

Aquí no hay método iterativo. El número de chapas y la separación entre paneles
son datos de entrada, tomados de la Tabla 9, y las dos velocidades salen
directamente de la ecuación (3), w = G/(rho*N*f_ch). La longitud necesaria sale
de la ec. (26) y las pérdidas de carga de la ec. (19), que se comparan con las
que reporta el artículo.

Dos avisos sobre la interpretación de los datos:

  - La columna n_pl de la Tabla 9 es el número de CHAPAS, no de canales. Cada
    panel se forma soldando dos chapas (de ahí que la ec. (13) obligue a que
    n_pl sea par), así que los canales internos son n_pl/2. Con n_pl las
    densidades que implica la Tabla 9 saldrían de 250 a 390 kg/m3, imposibles
    para hidrocarburos líquidos; con n_pl/2 salen de 506 a 770 kg/m3 y decrecen
    de forma monótona con la temperatura. Ver CHAPAS_POR_PANEL.

  - En este caso de estudio el fluido CALIENTE tiene menor caudal y va al canal
    INTERNO, y el crudo (frío) al externo. Es al revés que en el caso 1, y es
    además la configuración para la que se dedujeron las ecs. (22)-(28), en la
    que "subíndice 1 = caliente" y "subíndice 1 = canal I" coinciden.

Correlaciones tomadas de O. Arsenyeva et al.
"""

from dataclasses import dataclass, replace
import formulas as form
import asignacion as asig
import geometría as geom

CHAPAS_POR_PANEL = 2      # n_pl de la Tabla 9 son chapas; los canales internos son n_pl/CHAPAS_POR_PANEL
PROPIEDADES_DE_REF_24 = False   # ponlo a True cuando sustituyas rho, mu, k y cp por los valores reales

# =============================================================================
# Caso de estudio 2: condiciones de proceso (Tabla 7) y diseños óptimos (Tabla 9)
# =============================================================================

# Las seis variantes resultaron ser geometría PPHE2. n_pl y b son los de la
# Tabla 9; G, t_in y t_out los de la Tabla 7 (G venía en t/h, aquí en kg/s).
# La pérdida de carga admisible es de 100 kPa en los dos lados.
#
# PROPIEDADES: el artículo las toma de su referencia [24] y NO las publica. Las
# de aquí son provisionales: rho está despejada de la propia Tabla 9, de modo
# que la comparación de velocidades es circular mientras no metas las reales, y
# mu, k y cp son valores tipo de hidrocarburo, puestos solo para que el script
# corra. Sustitúyelos y pon PROPIEDADES_DE_REF_24 = True.

_MU, _K, _CP = 3.0e-4, 0.11, 2400.0      # provisionales, iguales para todo

VARIANTES = [
    dict(HE=1, n_pl=24, b=1.5e-3,
         int=dict(fluido="Producto 1", G=37.000/3.6, t_in=177.4, t_out=155.6, del_P=100e3,
                  rho=769.6, mu=_MU, k=_K, cp=_CP),
         ext=dict(fluido="Crudo",      G=107.339/3.6, t_in=150.0, t_out=157.0, del_P=100e3,
                  rho=720.4, mu=_MU, k=_K, cp=_CP),
         t9=dict(F=15.85, L_F=2.25, U=2863, w_int=1.943, w_ext=5.370, dP_int=1.000e5, dP_ext=5.417e4)),
    dict(HE=2, n_pl=24, b=1.0e-3,
         int=dict(fluido="Producto 2", G=39.200/3.6, t_in=228.3, t_out=145.4, del_P=100e3,
                  rho=710.8, mu=_MU, k=_K, cp=_CP),
         ext=dict(fluido="Crudo",      G=107.339/3.6, t_in=81.4, t_out=117.5, del_P=100e3,
                  rho=781.1, mu=_MU, k=_K, cp=_CP),
         t9=dict(F=9.43, L_F=1.35, U=2969, w_int=2.229, w_ext=6.271, dP_int=1.000e5, dP_ext=5.437e4)),
    dict(HE=3, n_pl=20, b=2.0e-3,
         int=dict(fluido="Residuo 1",  G=26.327/3.6, t_in=275.0, t_out=211.0, del_P=100e3,
                  rho=622.5, mu=_MU, k=_K, cp=_CP),
         ext=dict(fluido="Crudo",      G=107.339/3.6, t_in=161.2, t_out=178.5, del_P=100e3,
                  rho=655.4, mu=_MU, k=_K, cp=_CP),
         t9=dict(F=9.09, L_F=1.67, U=2125, w_int=2.051, w_ext=5.853, dP_int=1.000e5, dP_ext=3.440e4)),
    dict(HE=4, n_pl=58, b=0.0,
         int=dict(fluido="Producto 3", G=102.426/3.6, t_in=289.9, t_out=270.2, del_P=100e3,
                  rho=511.9, mu=_MU, k=_K, cp=_CP),
         ext=dict(fluido="Crudo",      G=82.422/3.6, t_in=224.6, t_out=250.5, del_P=100e3,
                  rho=575.8, mu=_MU, k=_K, cp=_CP),
         t9=dict(F=12.71, L_F=0.83, U=3316, w_int=3.346, w_ext=5.779, dP_int=1.000e5, dP_ext=2.677e4)),
    dict(HE=5, n_pl=62, b=0.5e-3,
         int=dict(fluido="Producto 4", G=102.426/3.6, t_in=289.9, t_out=264.4, del_P=100e3,
                  rho=548.8, mu=_MU, k=_K, cp=_CP),
         ext=dict(fluido="Crudo",      G=120.755/3.6, t_in=224.6, t_out=247.5, del_P=100e3,
                  rho=617.4, mu=_MU, k=_K, cp=_CP),
         t9=dict(F=19.13, L_F=1.10, U=2937, w_int=2.920, w_ext=4.708, dP_int=1.000e5, dP_ext=2.519e4)),
    dict(HE=6, n_pl=24, b=1.0e-3,
         int=dict(fluido="Residuo 2",  G=24.372/3.6, t_in=346.0, t_out=300.0, del_P=100e3,
                  rho=505.9, mu=_MU, k=_K, cp=_CP),
         ext=dict(fluido="Crudo",      G=82.422/3.6, t_in=269.3, t_out=284.0, del_P=100e3,
                  rho=548.9, mu=_MU, k=_K, cp=_CP),
         t9=dict(F=11.79, L_F=1.85, U=3050, w_int=1.947, w_ext=6.852, dP_int=1.000e5, dP_ext=3.978e4)),
]

# =============================================================================
# Estructuras de datos
# =============================================================================



@dataclass(frozen=True)
class PPHEResult:
    nombre: str
    N: float
    d_h_int: float
    u_m_int: float
    Pr_int: float
    Re_int: float
    Nu_int: float
    f_int: float
    h_int: float
    d_h_ext: float
    u_m_ext: float
    Pr_ext: float
    Re_ext: float
    Nu_ext: float
    f_ext: float
    h_ext: float
    U: float
    L_F: float
    F: float
    del_P_int: float
    del_P_ext: float
    t9: dict


def imprimir_resultados(resultados):
    """Da formato a la impresión por consola de los resultados."""
    encabezado = (
        f"{'Caso':<8}{'N':>6}{'Dh (in) [m]':>15}{'Dh (out) [m]':>15}{'u (in) [m/s]':>15}{'u (out) [m/s]':>15}{'Pr (in)':>10}{'Pr (out)':>10}{'Re (in)':>15}{'Re (out)':>15}{'Nu (in)':>12}{'Nu (out)':>12}{'f (in)':>12}{'f (out)':>12}{'h (in) [W/m2K]':>16}{'h (out) [W/m2K]':>16}{'U [W/m2K]':>14}{'L_F [m]':>12}{'dP (in) [kPa]':>15}{'dP (out) [kPa]':>15}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for res in resultados:
        print(
            f"{res.nombre:<8}{res.N:6.0f}{res.d_h_int:15.6e}{res.d_h_ext:15.6e}{res.u_m_int:15.4f}{res.u_m_ext:15.4f}{res.Pr_int:10.4f}{res.Pr_ext:10.4f}{res.Re_int:15.2e}{res.Re_ext:15.2e}{res.Nu_int:12.4f}{res.Nu_ext:12.4f}{res.f_int:12.4f}{res.f_ext:12.4f}{res.h_int:16.2f}{res.h_ext:16.2f}{res.U:14.2f}{res.L_F:12.4f}{res.del_P_int/1e3:15.2f}{res.del_P_ext/1e3:15.2f}"
        )


def imprimir_comparacion(resultados):
    """Contrasta lo calculado con lo que reporta la Tabla 9 del artículo."""
    encabezado = (
        f"{'Caso':<8}{'u_in':>9}{'T9':>9}{'u_out':>9}{'T9':>9}{'L_F [m]':>10}{'T9':>9}"
        f"{'U':>10}{'T9':>9}{'F [m2]':>10}{'T9':>9}{'dP_in [kPa]':>13}{'T9':>9}{'dP_out':>9}{'T9':>9}"
    )
    print('\n' + encabezado)
    print('-' * len(encabezado))
    for r in resultados:
        t = r.t9
        print(
            f"{r.nombre:<8}{r.u_m_int:9.3f}{t['w_int']:9.3f}{r.u_m_ext:9.3f}{t['w_ext']:9.3f}"
            f"{r.L_F:10.3f}{t['L_F']:9.3f}{r.U:10.0f}{t['U']:9.0f}{r.F:10.2f}{t['F']:9.2f}"
            f"{r.del_P_int/1e3:13.1f}{t['dP_int']/1e3:9.1f}{r.del_P_ext/1e3:9.1f}{t['dP_ext']/1e3:9.1f}"
        )



# =============================================================================
# Cálculos numéricos
# =============================================================================

if not PROPIEDADES_DE_REF_24:
    print("=" * 100)
    print("AVISO: propiedades de fluido PROVISIONALES. El artículo las toma de su ref. [24] y no las")
    print("publica. rho está despejada de la propia Tabla 9 (la comparación de velocidades es por tanto")
    print("circular) y mu, k y cp son valores genéricos de hidrocarburo. Todo lo que dependa de Nu, h,")
    print("U, L_F y las pérdidas de carga carece de valor hasta que metas los datos reales.")
    print("=" * 100)

resultados = []

for v in VARIANTES:
    nombre = f"HE{v['HE']}"
    props_fl_int, props_fl_ext = v["int"], v["ext"]

    # Todas las variantes de la Tabla 9 son geometría PPHE2; solo cambia la separación b
    i = replace(geom.PPHE2, b=v["b"])

    constantes_n = asig.asignacion_de_constantes(sT=i.s_t, s2L=i.s_2l, dsp=i.d_sp, h=i.b_i)

    """VELOCIDADES (ecuación (3): el caudal repartido entre N canales)"""

    N = v["n_pl"]/CHAPAS_POR_PANEL                                                                                              # Canales, no chapas
    f_ch_int = form.f_chI(i.b_i, i.w_pp, i.w_e)                                                                                 # Sección de paso (16)
    f_ch_ext = form.f_chE(i.b_i, i.b, i.w_pp, i.w_e)                                                                            # Sección de paso (17)
    u_m_int = form.velocidad_ec3(G=props_fl_int.get("G"), rho=props_fl_int.get("rho"), N=N, f_ch=f_ch_int)                       # Velocidad media
    u_m_ext = form.velocidad_ec3(G=props_fl_ext.get("G"), rho=props_fl_ext.get("rho"), N=N, f_ch=f_ch_ext)                       # Velocidad media

    """CÁLCULO INTERNO"""

    d_h_int = form.deI(i.b_i)                                                                                                   # Diámetro equivalente
    Re_int = form.Re(rho=props_fl_int.get("rho"), u=u_m_int, dh=d_h_int, mu=props_fl_int.get("mu"))                              # Número de Reynolds
    Pr_int = form.Pr(cp=props_fl_int.get("cp"), mu=props_fl_int.get("mu"), k=props_fl_int.get("k"))                              # Número de Prandtl
    Nu_int = form.NuI(n3=constantes_n[2], n4=constantes_n[3], n5=constantes_n[4], Re=Re_int, Pr=Pr_int)                          # Número de Nusselt
    f_int = form.fI(n1=constantes_n[0], Re=Re_int, n2=constantes_n[1])                                                           # Factor de fricción
    h_int = Nu_int * props_fl_int.get("k") / d_h_int                                                                             # Coef. de transferencia de calor


    """CÁLCULO EXTERNO"""

    d_h_ext = form.deE(i.b_i, i.b)                                                                                              # Diámetro equivalente
    Re_ext = form.Re(rho=props_fl_ext.get("rho"), u=u_m_ext, dh=d_h_ext, mu=props_fl_ext.get("mu"))                              # Número de Reynolds
    Pr_ext = form.Pr(cp=props_fl_ext.get("cp"), mu=props_fl_ext.get("mu"), k=props_fl_ext.get("k"))                              # Número de Prandtl
    f_ext = form.fE(Re=Re_ext)                                                                                                  # Factor de fricción
    Nu_ext = form.NuE(Re=Re_ext, Pr=Pr_ext, f=f_ext)                                                                            # Número de Nusselt
    h_ext = Nu_ext * props_fl_ext.get("k") / d_h_ext                                                                             # Coef. de transferencia de calor


    """LONGITUD Y PÉRDIDAS DE CARGA"""

    U_global = form.U(h_int, props_fl_int.get("k"), d_h_int, Nu_int,
                      h_ext, props_fl_ext.get("k"), d_h_ext, Nu_ext, i.delta_pp)                                                # Coef. global (29)
    NTU_0 = form.NTU(props_fl_int.get("t_in"), props_fl_int.get("t_out"),
                     props_fl_ext.get("t_in"), props_fl_ext.get("t_out"))                                                       # NTU mínimo (23)
    L_F = form.LF_thermal(i.b_i, NTU_0, props_fl_int.get("cp"), u_m_int,
                          props_fl_int.get("rho"), U_global, i.Fx)                                                              # Longitud necesaria (26)
    F = v["n_pl"] * L_F * (i.w_pp - 2*i.w_e) * i.Fx                                                                             # Área de intercambio (5, 25)
    del_P_int = form.perdida_de_carga(f_int, L_F, d_h_int, props_fl_int.get("rho"), u_m_int)                                     # Pérdida de carga (19)
    del_P_ext = form.perdida_de_carga(f_ext, L_F, d_h_ext, props_fl_ext.get("rho"), u_m_ext)                                     # Pérdida de carga (19)

    Q_int = props_fl_int.get("G")*props_fl_int.get("cp")*abs(props_fl_int.get("t_out")-props_fl_int.get("t_in"))
    Q_ext = props_fl_ext.get("G")*props_fl_ext.get("cp")*abs(props_fl_ext.get("t_out")-props_fl_ext.get("t_in"))
    print(f"\n--- {nombre}: {props_fl_int.get('fluido')} (interno) contra {props_fl_ext.get('fluido')} (externo) ---")
    print(f"n_pl={v['n_pl']} chapas -> N={N:.0f} canales,  b={v['b']*1e3:.1f} mm")
    print(f"u_int={u_m_int:.4f} m/s (Tabla 9: {v['t9']['w_int']:.3f}), u_ext={u_m_ext:.4f} m/s (Tabla 9: {v['t9']['w_ext']:.3f})")
    print(f"Balance Q_int/Q_ext = {Q_int/Q_ext:.4f}  <- solo será 1 con los cp reales de [24]")

    resultados.append(
        PPHEResult(
            nombre=nombre,
            N=N,
            d_h_int=d_h_int,
            u_m_int=u_m_int,
            Pr_int=Pr_int,
            Re_int=Re_int,
            Nu_int=Nu_int,
            f_int=f_int,
            h_int=h_int,
            d_h_ext=d_h_ext,
            u_m_ext=u_m_ext,
            Pr_ext=Pr_ext,
            Re_ext=Re_ext,
            Nu_ext=Nu_ext,
            f_ext=f_ext,
            h_ext=h_ext,
            U=U_global,
            L_F=L_F,
            F=F,
            del_P_int=del_P_int,
            del_P_ext=del_P_ext,
            t9=v["t9"],
        )
    )

imprimir_resultados(resultados)
imprimir_comparacion(resultados)
