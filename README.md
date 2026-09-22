# receptor_solar_PPHE

Diseño de un receptor solar basado en intercambiadores de tipo "pillow-plate"
para comparar sus prestaciones con otras alternativas de intercambiadores de
calor.

## Página interactiva

**https://pablogcguc3m.github.io/receptor_solar_PPHE/**

Reproduce el modelo nodal en el navegador y deja mover sus parámetros para ver
cómo cambia la distribución de calor: el mapa de flujo incidente, la
temperatura del aire y la de la pared sobre la placa, el perfil de salida en
`y = L` y el recorrido de cada columna. No necesita servidor ni instalar nada.

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

HTML, CSS y JavaScript sin dependencias. `fisica.js` es una traducción del
núcleo de cálculo de `modelo_receptor/`, validada nodo a nodo contra el
original. Se publica sola al empujar cambios a `main`.
