"""Este script recoge todas las fórmulas"""



from math import sqrt, log

# =============================================================================
# Números adimensionales
# =============================================================================

def fI(n1, Re, n2):
    """Factor de fricción para el CANAL INTERNO (ley de potencia)."""
    return n1 * Re ** n2

def fE(Re):
    """Factor de fricción para el CANAL EXTERNO (ley de potencia, se puede asumir cierta independientemente de la geometría con las constantes dadas)"""
    A = 2.187
    n = 0.356
    return A*Re**(-n)

def NuI(n3, n4, n5, Re, Pr):
    """Número de Nusselt para el CANAL INTERNO ley de potencia. Ecuación (6) del artículo de M.Piper"""
    return n3 * Re ** n4 * Pr ** n5

def NuE(Re, Pr, f):
    """Número de Nusselt para el CANAL EXTERNO (correlación experimental, se puede asumir cierta independientemente de la geometría)"""
    psi = 0.58
    return (psi*f/8*Re*Pr)/(1.07+12.7*(psi*f/8)**(1/2)*(Pr**(2/3)-1))

def Re(rho, u, dh, mu):
    """Número de Reynolds."""
    return (rho * u * dh) / mu


def Pr(cp, mu, k):
    """Número de Prandtl."""
    return (cp * mu) / k

def NTU(t1in,t1out,t2in,t2out):
    """Número de unidades de transferencia MÍNIMO (23) O.Arsenyeva.

    Referido al fluido 1 y en contracorriente. Se usan valores absolutos para
    que sea indiferente que el fluido 1 sea el caliente o el frío: los saltos
    terminales son |t1in-t2out| (extremo por el que entra el fluido 1) y
    |t1out-t2in| (extremo por el que sale).
    """
    dT_a = abs(t1in-t2out)
    dT_b = abs(t1out-t2in)
    if abs(dT_a-dT_b) < 1e-9*max(dT_a, dT_b):
        LMTD = (dT_a+dT_b)/2      # limite dT_a -> dT_b: log(1) = 0 y la expresion general se indetermina
    else:
        LMTD = (dT_a-dT_b)/log(dT_a/dT_b)
    return abs(t1in-t1out)/LMTD


# =============================================================================
# Velocidad y diámetros equivalentes
# =============================================================================


def deI(b_i):
    """Diámetro equivalente interno (14) O.Arsenyeva"""
    return 2 * (b_i / sqrt(2)) 

def deE(b_i, b):
    """Diámetro equivalente externo (15) O.Arsenyeva"""
    return 2 * ((b_i+b)-b_i/sqrt(2))

def f_chI(b_i, w_pp, w_e):
    """Sección de paso del canal INTERNO (16) O.Arsenyeva. Es deI/2 por la anchura útil."""
    return (b_i/sqrt(2))*(w_pp-2*w_e)

def f_chE(b_i, b, w_pp, w_e):
    """Sección de paso del canal EXTERNO (17) O.Arsenyeva. Es deE/2 por la anchura útil."""
    return ((b_i+b)-b_i/sqrt(2))*(w_pp-2*w_e)

def velocidad_ec3(G, rho, N, f_ch):
    """Velocidad media a partir del gasto másico (3) O.Arsenyeva.

    G es el gasto de toda la corriente y N el número de canales entre los que se
    reparte, de forma que N*f_ch es la sección de paso total de ese lado.
    """
    return G/(rho*N*f_ch)

def perdida_de_carga(f, L_F, d_e, rho, w, zeta_DZ=1.5):
    """Pérdida de carga en el canal (19) O.Arsenyeva, zonas de distribución incluidas.

    Es la inversa de LF_ploss: aquella despeja la longitud conocida la pérdida de
    carga, y ésta la pérdida de carga conocida la longitud.
    """
    return f*L_F/d_e*rho*w**2/2 + zeta_DZ*rho*w**2

def wI(del_P, rho, f, NTU, cp, w, U, Fx):
    """Cálculo de la velocidad interna según la fórmula (28) del artículo de O.Arsenyeva.

    OJO: el artículo imprime 8*sqrt(2) en el denominador, pero esa constante es
    incompatible con el resto de sus propias ecuaciones. Igualando las dos
    expresiones de L_F (LF_ploss y LF_thermal) y sustituyendo de1 = 2*b_i/sqrt(2)
    el b_i se cancela y el denominador queda 4*sqrt(2)*de1/b_i = 8. Con 8*sqrt(2)
    las dos longitudes difieren siempre en un factor sqrt(2); con 8 coinciden.
    """
    return sqrt(del_P/(rho)*1/(1.5+(f*NTU*cp*w*rho)/(8*U*Fx)))

def wE(wI, cp1, cp2, rho1, rho2, t1in, t1out, t2in, t2out, de1, de2, wpp, we):
    """Velocidad media por el conducto externo según la fórmula (30) de O. Arsenyeva."""
    fch1 = de1/2*(wpp-2*we)
    fch2 = de2/2*(wpp-2*we)
    return wI*(cp1*rho1*(t1in-t1out)*fch1)/(cp2*rho2*(t2out-t2in)*fch2)

def U(h1,k1,de1,Nu1,h2,k2,de2,Nu2,delta_w):
    """Coeficiente global de transmisión de calor (29) de O.Arsenyeva."""
    h1=Nu1*k1/de1
    h2=Nu2*k2/de2
    lambda_w=16  # W/(m*K)
    inv_U = 1/h1 + 1/h2 + delta_w/lambda_w
    return 1/inv_U

def LF_ploss(de1,f,del_P,rho,w):
    """Longitud para desarrollar la pérdida de carga"""
    return 2*de1/f*(del_P/(rho*w**2)-1.5)

def LF_thermal(b_i,NTU,cp,w,rho,U,Fx):
    """Longitud para desarrollar la transferencia de calor"""
    return (b_i*NTU*cp*w*rho)/(2*sqrt(2)*U*Fx)