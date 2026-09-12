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
    LMTD = ((t1in-t2out)-(t1out-t2in))/log((t1in-t2out)/(t1out-t2in))
    return (t1in-t1out)/LMTD


# =============================================================================
# Velocidad y diámetros equivalentes
# =============================================================================

def wI(del_P, rho, f, NTU, cp, w, U, Fx):
    """Cálculo de la velocidad interna según la fórmula (28) del artículo de O.Arsenyeva."""
    return sqrt(del_P/(rho)*1/(1.5+(f*NTU*cp*w*rho)/(8*sqrt(2)*U*Fx)))

def velocidad_media_ext(G, rho, b_i, b, w_pp, w_e):
    """Velocidad media por el conducto externo."""
    area_flujo = (b_i + b - b_i / sqrt(2)) * (w_pp - 2 * w_e)
    return (G / (rho * area_flujo))  # Se aplica un factor de corrección


def diametro_hidraulico_interno(b_i):
    """Diámetro hidráulico interno."""
    return 2 * (b_i / sqrt(2))  # Se aplica un factor de corrección

def diametro_hidraulico_externo(b_i, b):
    """Diámetro hidráulico externo."""
    return 2 * ((b_i+b)-b_i/sqrt(2))
