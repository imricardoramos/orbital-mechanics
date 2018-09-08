# MODELAMIENTO DE DINÁMICA ORBITAL DE CUBESAT 3U PARA DETERMINACIÓN DE COSTOS PROPULSIVOS, ENERGÉTICOS Y TEMPORALES EN MANIOBRAS ORBITALES DE BAJO EMPUJE PREDETERMINADOS
### Explicación del código implementado a la fecha:
Se han implementado 5 archivos:
- <strong>constants.py:<strong> en éste código se definen todas las constantes a ser usadas en el resto del código como las constantes gravitacionales estándares de la Tierra, el Sol y la Luna ($\mu_E$, $\mu_M$, $\mu_S$), el radio de la Tierra $R_E$, la constante de gravedad $g_0$ y las masas de la Tierra, la Luna y el Sol ($M_E$, $M_M$, $M_S$), entre otras.
- <strong>coordinates.py:</strong> en éste código se definen funciones que son de ayuda para transformar las coordenadas de keplerianas a cartesianas, de cartesianas a keplerianas a cartesianas, de keplerianas a equinocciales, de equinocciales a keplerianas, de cartesianas a queinocciales y de equinocciales a cartesianas.
- <strong>models.py:</strong> en éste código se definen clases y funciones que conforman los modelos físicos a utilizar, como por ejemplo el modelo de densidad atmosférica, el modelo de posición Lunar y solar, y finalemente los modelos utilizados en el satélite, como los modelos de propulsores, los modelos de paneles solares, baterías, energía y los modelos de satélite en sí (incorporando cosas como masa, coeficiente balístico, etc). La idea es que acá se vayan agregando los modelos que sean distintos o más complejos en el tiempo.
- <strong>maneuvers.py:</strong> en éste código se encuentra las funciones que realizan el trabajo de calcular la dinámica orbital, y las funciones con las que interactúa el usuario para poder añadir las perturbaciones que desee y propagar el movimiento por el tiempo que se desee. En éste caso también se definen varios métodos de integración del movimiento del satélite, como el método de Conwell, el de Gauss y el de Betts. Un objeto \textit{Maneuver}, va calculando la trayectoria aplicada sobre él y guardandola en un historial que sirve para ser utilizado más tarde.
- <strong>helpers.py:</strong> en éste código se definen funciones auxiliares o misceláneas que tienen como propósito facilitar la interacción del usuario. En este momento sólo se tiene un extractor de elementos a partir de un TLE un una funcíón para facilitar la mantención del aspect ratio al graficar en 3D.

### Modo de uso:
En primer lugar se deben definir el modelo de <i>Spacecraft</i>, los elementos orbitales asociados a él y la fecha a la cual pertenecen dichos elementos orbitales.
Si bien en el ejemplo de abajo se obtienen a partir de un TLE, se pueden definir manualmente.
`coe` corresponde a los elementos orbitales expresados de manera kepleriana y en una lista de orden `[a,e,i,omega,Omega,M]`. La fecha simplemente es un `datetime`.
Éstos parametos son pasados al constructor del objeto Maneuvers, para definir el estado inicial de la maniobra.
```python
coe,date = helpers.parseTle("suchai0.tle")
satellite = Cubesat()
maneuvers = Maneuvers(coe,satellite,date)
```
Luego podemos agregar las perturbaciones al objeto maneuvers para agregar las perturbaciones (hasta ahora se ha implementado `atmosphere`, `solar_pressure`, `moon_gravity`, `sun_gravity`, `J2` y ``thrust`)
El método `propagate2` propaga en el tiempo y acepta el tiempo en segundos, con un mínimo de 60 segundos.
```python
# Add solar pressure and atmospheric drag perturbations to maneuver
maneuvers.addPerturbation("solar_pressure")
maneuvers.addPerturbation("atmosphere")
# Propagate 1 day 
maneuvers.propagate2(60*60*24*1)
# Start thrust
maneuvers.addPerturbation("thrust")
# Propagate for 18 days
maneuvers.propagate2(60*60*18)
# Stop thrust
maneuvers.removePerturbation("thrust")
# Propagate for 1 day
maneuvers.propagate2(60*60*24*1)
```
