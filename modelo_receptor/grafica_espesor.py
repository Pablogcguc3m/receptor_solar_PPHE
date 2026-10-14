"""Gráfica de t = p * sqrt(P / (3 S)) frente a P, para límite elástico y resistencia última.

    t = p * sqrt( P / (3 * S) )

P es la variable (eje x, en bar). S toma dos valores: el del límite elástico
(2.9e5 kPa) y el de la resistencia última (5.8e5 kPa). p = 116.108 mm es constante.
El texto se compone con LaTeX (Latin Modern), así que hace falta pdflatex instalado.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

p = 116.108  # [mm]

plt.rcParams.update({"text.usetex": True,  # Latin Modern vía LaTeX (requiere pdflatex)
                     "font.family": "serif",
                     "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath,bm}"})
S_ELASTICO = 2.9e5
S_ULTIMA = 5.8e5

KPA_POR_BAR = 100.0  # 1 bar = 100 kPa
P_MIN = 0.0  # rango de la variable P [bar]
T_MAX = 5.0    # límite superior del eje de ordenadas [mm]
P_MAX = 30.0  # [bar]
P_MARCA_MPA = 0.9  # presión marcada con la línea roja vertical [MPa]


def tex_sci(x):
    m, e = f"{x:.1e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def espesor(P, S, p=p):
    return p * np.sqrt(P / (3.0 * S))


SALIDA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "informes", "figuras", "espesor_vs_P.png")


def main(salida=SALIDA):
    P = np.linspace(P_MIN, P_MAX, 400)  # [bar]
    P_kPa = P * KPA_POR_BAR  # S está en kPa
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(P, espesor(P_kPa, S_ELASTICO), label=f"Límite elástico ($S = {tex_sci(S_ELASTICO)}$ kPa)")
    ax.plot(P, espesor(P_kPa, S_ULTIMA), label=f"Resistencia última ($S = {tex_sci(S_ULTIMA)}$ kPa)")
    ax.set_xlabel(r"\textbf{$\bm{P}$ [bar]}", labelpad=12)
    ax.set_ylabel(r"\textbf{$\bm{t = p\,\sqrt{P/(3S)}}$ [mm]}", labelpad=12)
    ax.axvline(P_MARCA_MPA * 10.0, color="red", linestyle="--", linewidth=0.8,
               label=f"Pérdida de carga ($\\Delta P = {P_MARCA_MPA * 10:g}$ bar)")  # 1 MPa = 10 bar
    ax.set_xlim(P_MIN, P_MAX)
    ax.set_ylim(0, T_MAX)
    ax.set_title(f"$t$ en función de $P$  ($p = {p}$ mm)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    if salida:
        os.makedirs(os.path.dirname(salida), exist_ok=True)
        fig.savefig(salida, dpi=200)
    plt.show()


if __name__ == "__main__":
    main()
