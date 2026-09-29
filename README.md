# receptor_solar_PPHE

Diseño de un receptor solar basado en intercambiadores de tipo "pillow-plate"
para comparar sus prestaciones con otras alternativas de intercambiadores de
calor.

## Páginas interactivas

Todas corren en el navegador, sin servidor ni instalar nada.

- **https://pablogcguc3m.github.io/receptor_solar_PPHE/** — la portada:
  animaciones ilustrativas de un receptor de pillow plates (la placa que se
  calienta bajo una mancha que sigue al puntero, el corte de un pillow plate
  inflándose y la planta del octógono con el campo de helióstatos). No calcula
  nada: enlaza las dos páginas siguientes.
- **https://pablogcguc3m.github.io/receptor_solar_PPHE/aire/** — la placa con
  aire (`receptor.py`), para demostraciones: el mapa de flujo incidente, la
  temperatura del aire y la de la pared, el perfil de salida en `y = L` y el
  recorrido de cada columna.
- **https://pablogcguc3m.github.io/receptor_solar_PPHE/octogonal/** — el caso
  de estudio final (`receptor_octogonal.py`): el receptor octogonal de sal
  solar, con la geometría del panel dentro de las correlaciones de Piper, las
  pérdidas y la mancha como parámetros. Da el gasto, los campos de sal, pared
  expuesta y pared mojada, el recorrido de la sal y el reparto placa a placa.

## Estructura

### `modelo_receptor/` — el receptor solar

- **`props_fluido.py`**: propiedades del fluido de trabajo en función de la
  temperatura. Aire como gas ideal, con `cp`, `k` y `mu` ajustados a CoolProp
  por polinomios de grado 5 entre 250 y 1600 K. Incluye la entalpía integrada.
- **`mapa.py`**: flujo solar incidente `q''(x, y)` sobre la placa, repartido
  como una gaussiana. La potencia se integra sobre la placa con la función
  error, de modo que el reparto por nodos conserva la energía exactamente.
- **`receptor.py`**: balance de energía porción a porción, con una porción por
  celda del patrón de soldaduras. Pérdidas por convección y radiación,
  conducción por la chapa y convección interna con el coeficiente del canal de
  pillow plate, descontando el área de las soldaduras. La mezcla entre columnas
  se regula con el factor `lambda`. Da el campo de temperaturas, el perfil
  `T(x)` en `y = L`, la pérdida de carga y una tabla de rendimiento frente a la
  temperatura de entrada.
- **`receptor_octogonal.py`**: el caso de estudio final. Receptor octogonal de
  8 pillow plates de 10.5 m con sal solar (propiedades de Zavoico), resuelto
  por simetría sobre medio receptor: 4 placas en serie, con mezcla completa
  entre placas. Fija la entrada (290 C) y la salida (565 C) y despeja el gasto.
  Mismo balance por porción que `receptor.py`, vectorizado por filas. Da los
  campos de temperatura de cada placa, la pérdida de carga y, opcionalmente,
  un barrido de `lambda` (`BARRIDO_LAMBDA`). `receptor.py` se mantiene como
  modelo de demostración con aire.
- **`documentacion_receptor.tex`**: documentación del módulo anterior.

### `modelo_termohidraulico/` — el canal de la pillow plate

- **`asignacion.py`**: asigna las constantes de las leyes de potencia según el
  caso de geometría que se tenga.
- **`formulas.py`**: agrupa todas las fórmulas necesarias para calcular los
  parámetros termohidráulicos.
- **`geometria.py`**: las geometrías de panel a estudiar (PPHE1, PPHE2, PPHE3).
- **`iteracion_vel.py`**: cálculo iterativo de la velocidad del fluido, que
  viene dada por una ecuación implícita.
- **`modelo_f3.py`**: calcula los resultados finales.

### `web/` — la página

HTML, CSS y JavaScript sin dependencias. `index.html` es la portada,
`aire/` y `octogonal/` las dos páginas de modelo. `aire/fisica.js` es una traducción del
núcleo de cálculo de `receptor.py` y `octogonal/fisica_octogonal.js`, la de
`receptor_octogonal.py`; las dos validadas nodo a nodo contra el original. Se
publica sola al empujar cambios a `main`.
