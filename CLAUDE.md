# CLAUDE.md

Guía para trabajar en este repositorio. Es el código del TFG de un receptor
solar de torre hecho con intercambiadores pillow-plate (PPHE). El README da la
descripción de cada fichero; aquí va lo que no se ve leyendo un solo fichero.

## Idioma y estilo

- Todo en **español**: nombres de variables, docstrings, comentarios, salidas
  por consola, mensajes de commit e informes.
- Los docstrings de módulo son largos y argumentados, con apartados en
  mayúsculas (`GEOMETRÍA.`, `POR QUÉ NO SUTHERLAND.`, `HIPÓTESIS.`...). Cada
  constante lleva su unidad entre corchetes y, si su valor es discutible, una
  nota debajo que lo justifica con números. Mantener ese estilo: el código se
  defiende en la memoria del TFG, y la justificación forma parte del código.
- Cada ecuación cita su número en el artículo de origen, por defecto
  O. Arsenyeva et al., Appl. Therm. Eng. 147 (2019) 579-591: `(19)`, `ec. (3)`.
- Los `.py` llevan tildes. Los `.js`, el `index.html` de cada página y el
  workflow de GitHub se escriben **sin tildes** (ASCII) en los comentarios.
- Commits: una frase corta en español, en presente o impersonal, terminada en
  punto (`Añade sección...`, `Cambio correlación de nusselt`).

## Cómo se ejecuta

No hay paquete ni instalación. Los imports son planos (`import mapa`,
`import formulas`), así que **cada script se lanza desde su propia carpeta**:

```powershell
cd modelo_termohidraulico; python modelo_f3.py      # contraste con las tablas 5 y 9 de Arsenyeva
cd modelo_receptor; python receptor.py              # placa de aire 1.5 x 1.5 m (~3 s)
cd modelo_receptor; python receptor_octogonal.py    # caso final, sal solar (~5 s)
cd modelo_receptor; python receptor_octogonal_sabourishirazi.py
cd modelo_receptor; python props_fluido.py          # tabla y validación frente a CoolProp
cd modelo_receptor; python mapa.py                  # comprobaciones del mapa gaussiano
```

- Dependencias: Python 3.12, numpy y matplotlib. CoolProp solo hace falta para
  `props_fluido.validar_con_coolprop()`; el modelo corre sin él.
- Los scripts abren figuras con `plt.show()`. Para ejecutarlos sin ventanas,
  `MPLBACKEND=Agg`.
- `grafica_espesor.py` usa `text.usetex` (necesita pdflatex) y escribe en
  `informes/figuras/`.
- **No hay tests.** La verificación son los bloques `if __name__ == "__main__"`:
  contrastes con las tablas publicadas, con CoolProp y balances que deben
  cerrar. Tras un cambio de física, ejecutar el script afectado y comparar sus
  números con los de antes.

### La web

`web/` es HTML + JS sin dependencias ni compilación. Para verla en local está
la configuración `web` de `.claude/launch.json` (`python -m http.server 8765
--directory web`). Los dos núcleos de cálculo exportan con `module.exports`
(`RCP` en `aire/fisica.js`, `ROC` en `octogonal/fisica_octogonal.js`), así que
se pueden cargar en Node para validarlos contra Python.

Al empujar a `main` con cambios en `web/`, `.github/workflows/pages.yml` la
publica en GitHub Pages. Ese workflow añade `?v=<commit>` al script con `sed`
buscando el literal exacto `src="fisica.js"` / `src="fisica_octogonal.js"`, y
falla si no lo encuentra: no cambiar esa etiqueta `<script>`.

## Arquitectura

```
modelo_termohidraulico/          canal de la pillow plate, propiedades constantes
  geometria.py   PPHE1/2/3 (Tabla 1), dataclass congelada PPHEGeometry
  asignacion.py  n1..n5 de las leyes de potencia (tablas de Piper) según la geometría
  formulas.py    ecuaciones de Arsenyeva (Re, Pr, Nu, f, de, f_ch, U, L_F, dp)
  modelo_f3.py   vía directa: velocidad del caudal, ec. (3)
  iteracion_vel.py  vía inversa: velocidad de la ec. implícita (28); nadie la llama ahora

modelo_receptor/                 receptor solar, propiedades función de T
  props_fluido.py  AIRE (gas ideal, polinomios ajustados a CoolProp) y SAL_SOLAR (Zavoico)
  mapa.py          flujo incidente gaussiano, integrado con erf por celda
  receptor.py      placa plana de aire: balance porción a porción, bisección anidada
  receptor_octogonal.py   caso final: 8 placas de sal, medio receptor, Newton vectorizado con numpy
  receptor_octogonal_sabourishirazi.py  subclase que solo cambia el Nusselt interno
```

- `receptor.py` carga `formulas`, `asignacion` y `geometria` del otro paquete
  **por ruta** con `_importar()`, porque `modelo_receptor/formulas.py` existe
  (vacío) y un `import formulas` cogería el de la carpeta equivocada. No
  sustituirlo por un `sys.path.append`.
- `receptor_octogonal.py` reutiliza de `receptor.py` las constantes, las
  pérdidas y la cadena del canal interno; `receptor_octogonal_sabourishirazi.py`
  hereda de `ReceptorOctogonal`. Un cambio en `receptor.py` se propaga a los
  dos.
- **Python es la referencia y la web es una traducción.** `web/aire/fisica.js`
  porta `receptor.py` (con `mapa`, `props_fluido` y el paquete termohidráulico)
  y `web/octogonal/fisica_octogonal.js` porta `receptor_octogonal.py` y la
  variante de Sabourishirazi (`cfg.nusselt`). Cualquier cambio de ecuaciones,
  constantes o resolvedor en Python hay que llevarlo al JS correspondiente y
  validarlo nodo a nodo contra el original.

## Convenios del modelo que es fácil romper

- **Unidades SI y temperaturas en kelvin** dentro de `props_fluido` y de los
  modelos del receptor. Las condiciones de proceso se escriben en Celsius
  sumando `pf.CERO_CELSIUS` (alias `C` en `receptor_octogonal.py`).
  `modelo_termohidraulico` trabaja en Celsius con propiedades constantes.
- **Subíndices**: `I`/`E` son canal interno/externo; `1`/`2` son los del
  artículo (1 caliente, 2 frío). Solo coinciden si el caliente va por dentro.
- **Desviaciones deliberadas del artículo de Arsenyeva**, justificadas en los
  docstrings de `formulas.py`: `wI` usa 8 donde el artículo imprime 8·√2, y
  `NTU` toma valores absolutos. No "corregirlas" de vuelta.
- Los coeficientes que da `asignacion.py` (tablas de Piper) no coinciden con los
  que Arsenyeva tabula para PPHE1-3 (hasta 51 % en ζ). Es conocido y está
  documentado; `asignacion` lanza `ValueError` fuera de los rangos geométricos.
- En las tablas del artículo `n_pl` es el número de **chapas**; los canales
  internos son `n_pl/2` (`CHAPAS_POR_PANEL`).
- Las geometrías son `frozen=True`: se modifican con `dataclasses.replace()`,
  igual que el `Receptor` para barridos.
- La malla del receptor la fija el panel: una porción por celda de soldadura,
  `dx = s_T`, `dy = s_2l`, sin soldaduras de borde (`w_e = 0`). El flujo se
  pierde por toda el área `A`, pero al fluido solo pasa por `A_f` (sin puntos
  de soldadura).
- La mezcla entre columnas `LAMBDA` (1 adiabáticas, 0 mezcla completa) pondera
  **entalpías**, no temperaturas, para conservar la energía.
- **Fuera de rango**: las correlaciones de Piper valen para Re 1000-8000
  (`RANGO_RE_PIPER`) y la sal circula a Re 4.5e4-1.4e5; por defecto se
  extrapola y `comparacion_rango()` lo cuantifica. `props_fluido` lanza
  `ValueError` fuera del rango de ajuste (sal 260-600 C, aire 250-1600 K)
  salvo con `EXTRAPOLAR = True`.
- Los interruptores del caso de estudio son constantes al principio de
  `receptor_octogonal.py`: `THETA_ESTRELLA` (0 octógono, hasta 67.5° estrella),
  `BARRIDO_LAMBDA`, `COMPARAR_RANGO`, `DIAMETRO_CIRCUNSCRITO`.

## Informes

Al terminar cada cambio de código se escribe un informe LaTeX en `informes/`
con nombre `informe_AAAA-MM-DD_tema.tex` (article, babel spanish, booktabs,
compilable con pdflatex): qué se cambió, por qué, cómo y con qué resultados,
con las tablas numéricas. `informes/` y `.claude/` están en `.gitignore` y no
se suben. `modelo_receptor/documentacion_receptor.tex` sí se versiona, pero no
su PDF.

## Cabos sueltos conocidos

- El docstring de `modelo_f3.py` habla de `modelo_iterativo.py`, que no existe:
  la vía iterativa es `iteracion_vel.py`.
- `modelo_receptor/formulas.py` está vacío. Mientras exista, `_importar()` es
  imprescindible.
