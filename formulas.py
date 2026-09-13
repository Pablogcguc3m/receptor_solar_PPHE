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
    """Número de unidades de transferencia MÍNIMO (23) O.Arsenyeva"""
    LMTD = ((t1in-t2out)-(t1out-t2in))/log((t1in-t2out)/(t1out-t2in))
    return (t1in-t1out)/LMTD


# =============================================================================
# Velocidad y diámetros equivalentes
# =============================================================================


def deI(b_i):
    """Diámetro equivalente interno (14) O.Arsenyeva"""
    return 2 * (b_i / sqrt(2)) 

def deE(b_i, b):
    """Diámetro equivalente externo (15) O.Arsenyeva"""
    return 2 * ((b_i+b)-b_i/sqrt(2))

def wI(del_P, rho, f, NTU, cp, w, U, Fx):
    """Cálculo de la velocidad interna según la fórmula (28) del artículo de O.Arsenyeva."""
    return sqrt(del_P/(rho)*1/(1.5+(f*NTU*cp*w*rho)/(8*sqrt(2)*U*Fx)))

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