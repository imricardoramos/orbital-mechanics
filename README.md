# Modelamiento de dinámica orbital de cubesat 3U para determinación de costos propuslivos, energéticos y temporales en maniobras orbitales de bajo empuje predeterminados
### Explicación del código implementado a la fecha:
Se han implementado 5 archivos:
- <strong>constants.py:</strong> en éste código se definen todas las constantes a ser usadas en el resto del código como las constantes gravitacionales estándares de la Tierra, el Sol y la Luna, el radio de la Tierra, la constante de gravedad y las masas de la Tierra, la Luna y el Sol, entre otras.
- <strong>coordinates.py:</strong> en éste código se definen funciones que son de ayuda para transformar las coordenadas de keplerianas a cartesianas, de cartesianas a keplerianas a cartesianas, de keplerianas a equinocciales, de equinocciales a keplerianas, de cartesianas a queinocciales y de equinocciales a cartesianas.
- <strong>models.py:</strong> en éste código se definen clases y funciones que conforman los modelos físicos a utilizar, como por ejemplo el modelo de densidad atmosférica, el modelo de posición Lunar y solar, y finalemente los modelos utilizados en el satélite, como los modelos de propulsores, los modelos de paneles solares, baterías, energía y los modelos de satélite en sí (incorporando cosas como masa, coeficiente balístico, etc). La idea es que acá se vayan agregando los modelos que sean distintos o más complejos en el tiempo.
- <strong>maneuvers.py:</strong> en éste código se encuentra las funciones que realizan el trabajo de calcular la dinámica orbital, y las funciones con las que interactúa el usuario para poder añadir las perturbaciones que desee y propagar el movimiento por el tiempo que se desee. En éste caso también se definen varios métodos de integración del movimiento del satélite, como el método de Conwell, el de Gauss y el de Betts. Un objeto \textit{Maneuver}, va calculando la trayectoria aplicada sobre él y guardandola en un historial que sirve para ser utilizado más tarde.
- <strong>auxiliary.py:</strong> en éste código se definen funciones auxiliares o misceláneas que tienen como propósito facilitar la interacción del usuario. En este momento sólo se tiene un extractor de elementos a partir de un TLE un una funcíón para facilitar la mantención del aspect ratio al graficar en 3D.

### Modo de uso:
En primer lugar se debe declarar el modelo de <i>spacecraft</i>, los elementos orbitales asociados a él y la fecha a la cual pertenecen dichos elementos orbitales.
Si bien en el ejemplo de abajo se obtienen a partir de un TLE, se pueden definir manualmente.
`coe` corresponde a los elementos orbitales expresados de manera kepleriana y en una lista de orden `[a,e,i,omega,Omega,M]`. La fecha simplemente es un `datetime`.
Éstos parametos son pasados al constructor del objeto Maneuvers, para definir el estado inicial de la maniobra.
```python
coe,date = helpers.parseTle("suchai0.tle")
satellite = Cubesat("3U")
maneuvers = Maneuvers(coe,satellite,date)
```
Luego podemos agregar las perturbaciones al objeto maneuvers para agregar las perturbaciones (hasta ahora se ha implementado `atmosphere`, `solar_pressure`, `moon_gravity`, `sun_gravity`, `J2` y `thrust`)
El método `propagate` propaga en el tiempo y acepta el tiempo en segundos, y un timestep en segundos.
```python
# Add solar pressure and atmospheric drag perturbations to maneuver
maneuvers.addPerturbation("solar_pressure")
maneuvers.addPerturbation("atmosphere")
# Propagate 1 day 
maneuvers.propagate(60*60*24*1)
# Start thrust
maneuvers.addPerturbation("thrust")
# Propagate for 18 days
maneuvers.propagate2(60*60*24*18)
# Stop thrust
maneuvers.removePerturbation("thrust")
# Propagate for 1 day
maneuvers.propagate2(60*60*24*1)
```
En cada propagación se van guardando los datos en el historial de la maniobra:
```python
# Spacecraft distance vectors from Earth Center
maneuvers.history.r
# Spacecraft velocity vectors history
maneuvers.history.v
# Classical Orbital Elements History
maneuvers.history.coe
# Modified Equinoctial Elements History
maneuvers.history.mee
# Propellant Mass History
maneuvers.history.propMass
# Time elapsed history
maneuvers.history.t
# Datetime history
maneuvers.history.datetime
# Delta-V history
maneuvers.history.dv
```
Ejemplos de uso y resultados:
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Ascent.ipynb">Demo - Ascent</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Deorbiting.ipynb">Demo - Deorbiting</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Inclination Change.ipynb">Demo - Inclination Change</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Perturbations.ipynb">Demo - Perturbations</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Energy.ipynb">Demo - Energy</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Energy.ipynb">Demo - Energy</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Relative Motion.ipynb">Demo - Relative Motion</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Single Orbital Parameter Modification.ipynb">Demo - Single Orbital Parameter Modification</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Validation - Perturbations">Validation - Perturbations</a>
- <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Validation - STK">Validation - STK</a>
