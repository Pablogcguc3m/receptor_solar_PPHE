# receptor_solar_PPHE
Diseño de un receptor solar basado en intercambiadores de tipo "pillow-plate" para comparar sus prestaciones con otras alternativas de intercambiadores de calor.

La estructura de los scripts es la siguiente:
-asignación.py: Asigna valores de ley de potencia en según el caso de geometría que se tenga.
-formulas.py: Agrupa todas las fórmulas necesarias para calcular todos los parámetros termohidráulicos.
-modelo_termohidraulico.py: Calcula los resultados finales.
-iteracion_vel.py: Realiza un cálculo iterativopara calcular la velocidad del fluido, ya que esta viene dada por una ecuación implícita.  