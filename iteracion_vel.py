"""Este script calcula, mediante un método iterativo la velocidad en el canal interno y externo"""
from math import sqrt, log
from scipy.optimize import brentq
from formulas import (fI, fE, NuI, NuE, Re, Pr, NTU, deI, deE, wI, wE, U, LF_ploss, LF_thermal)


def resolver_velocidad(rho1, mu1, cp1, k1, n1_fric, n2_fric, n3_nu, n4_nu, n5_nu,
                        rho2, mu2, cp2, k2, t1in, t1out, t2in, t2out, del_P1, bi, b, wpp, we, Fx,
                        delta_w, w1_inicial=1.0, tol=1e-6, max_iter=200, relajacion=1.0):

    de1 = deI(bi)
    de2 = deE(bi, b)
    NTU0 = NTU(t1in, t1out, t2in, t2out)

    w1 = w1_inicial
    for it in range(max_iter):
        # --- canal I (con w1 actual) ---
        Re1 = Re(rho1, w1, de1, mu1)
        f1 = fI(n1_fric, Re1, n2_fric)
        Pr1 = Pr(cp1, mu1, k1)
        Nu1 = NuI(n3_nu, n4_nu, n5_nu, Re1, Pr1)
        h1 = Nu1 * k1 / de1

        # --- canal E (w2 depende de w1 vía Ec.30) ---
        w2 = wE(w1, cp1, cp2, rho1, rho2, t1in, t1out, t2in, t2out, de1, de2, wpp, we)
        Re2 = Re(rho2, w2, de2, mu2)
        f2 = fE(Re2)
        Pr2 = Pr(cp2, mu2, k2)
        Nu2 = NuE(Re2, Pr2, f2)
        h2 = Nu2 * k2 / de2

        U_actual = U(h1, k1, de1, Nu1, h2, k2, de2, Nu2, delta_w)

        w1_nuevo = wI(del_P1, rho1, f1, NTU0, cp1, w1, U_actual, Fx)
        w1_relajado = w1 + relajacion*(w1_nuevo - w1)

        if abs(w1_relajado - w1) < tol:
            w1 = w1_relajado
            print(f"Convergencia en {it+1} iteraciones. w1 = {w1:.5f} m/s")
            break
        w1 = w1_relajado
    else:
        print(f"AVISO: no convergió en {max_iter} iteraciones (último w1={w1:.5f})")

    # --- resultados finales, con w1 ya convergida ---
    LF_1 = LF_ploss(de1, f1, del_P1, rho1, w1)
    LF_2 = LF_thermal(bi, NTU0, cp1, w1, U_actual, Fx)
    print(f"w1={w1:.4f} m/s, w2={w2:.4f} m/s, Re1={Re1:.1f}, Re2={Re2:.1f}")
    print(f"LF (por presión)={LF_1:.4f} m, LF (por calor)={LF_2:.4f} m  <- deben coincidir")
    return w1, w2, LF_1, U_actual