# Modelling of the orbital mechanics of a 3U cubesat for the determination of propulsive, energetic and temporary costs in predetermined low thrust orbital maneuvers.
If you want, you can checkout <a href="https://ricardoramos.me/orbital-mechanics-simulator">the blog entry I wrote fot this project.<a>
    
### Explanation of the code written so far:
There are 5 files:
- <strong>constants.py:</strong> contains all the constants used in the rest of the code, like the standard gravitational constants of the Earth, the Moon and the Sun, the radius of the Earth, the gravitational constant and the masses of the Earth, the Moon and the Sun, among others.
- <strong>coordinates.py:</strong> contains functions which help in the transfomation of coordinates, from keplerian to cartesian, cartesian to equinoctial and equinoctial to cartesian.
- <strong>models.py:</strong> contains classes and functions which represent the physical models for use. For example the atmospheric density model, the lunar and solar positioning and finally the models used in the satellite units, like the thrusters, the solar panels, batteries, energy balance and even the satellite itself (incorporating stuff like mass, ballistic coefficient, etc). The idea is to keep adding models which are diverse and more complex over time.
- <strong>maneuvers.py:</strong> this is the core of the simulator and contains the functions which do the job of calculating the actual orbital dynamics, and the functions which the user interacts with to be able to add and remove perturbations, and propagate movement over time. In this case various methods of integration for the movement of the satellite have been defined, like the Cowell's method, Gauss with equinoctial coordinates (which I call Betts). The `Maneuver` object calculates the trajectory and logs it in a history class which can be further used.
- <strong>auxiliary.py:</strong> contains auxiliary or miscellaneous functions which don't classify under the other files. At the moment it only has an orbital elements extractor from a TLE, and a function to aid in keeping the aspect ratio of the 3D plot.

### Usage:
You must first declare the spacecraft model you wish to use, the initial orbital elements and the date they correspond to.
The example below defines the from a TLE, but they can be manually defined.
`coe` are the orbital elements expressed in keplerian form and in a list of the form `[a,e,i,omega,Omega,nu]`. 
The date is simply a `datetime` object.

The `Spacecraft` model internally defines a configuration for the thrusters, solar panels and batteries.
If you choose to use the thrusters, solar panels or batteries, you have to define the n in `Spacecraft.thruster`, `Spacecraft.solarPanles` and `Spacecraft.battery`. These can all be defined generically with `models.Thruster`, `models.solarPanels` and `models.Battery`, or any available model in <i>models.py</i>.

All these parameters are passed to the constructor of the `Maneuvers` object to define the initial state of the maneuver:

```python
coe,date = auxiliary.parseTle("suchai0.tle")
# Spacecraft Definition
satellite = models.Cubesat("3U")

# Thruster Definition (1mN thrust, 720s Isp)
satellite.thruster = models.Thruster(thrust=1e-3,isp=720)

# Solar Panels Definition
satellite.solarPanels = models.solarPanel()
satellite.solarPanels.area = 30e-2*30e-2
satellite.solarPanels.efficiency = 0.4

# Battery Definition
satellite.battery = models.NanoPowerBP4("2S-2P")

# Define maneuvers object
maneuvers = Maneuvers(coe,satellite,date)
```
We can then add the perturbations to the maneuvers object (so far only `atmosphere`, `solar_pressure`, `moon_gravity`, `sun_gravity`, `J2` and `thrust` have been implemented).
The `propagate` method propagates the time and accepts time in seconds, with a timestep in seconds.
```python
# Add solar pressure and atmospheric drag perturbations to maneuver
maneuvers.addPerturbation("solar_pressure")
maneuvers.addPerturbation("atmosphere")
# Propagate 1 day 
maneuvers.propagate(60*60*24*1, 60)
# Start thrust
maneuvers.addPerturbation("thrust")
# Propagate for 18 days
maneuvers.propagate(60*60*24*18, 60)
# Stop thrust
maneuvers.removePerturbation("thrust")
# Propagate for 1 day
maneuvers.propagate(60*60*24*1, 60)
```
To orient thrust, it is possible to define callback functions which accept `coe` (Classical Orbital Elements) as input, and output angles `alpha` and `beta` which define the orientation in the RSW reference system (also called RCN) of the satellite:
<img src="rswFrame.png"/>
```python
def alphaCallback(coe):
    # Alpha is such that is always pointing in the direction of velocity
    e = coe[1]
    nu = coe[5]
    alpha = np.arctan2(e*np.sin(nu),1+e*np.cos(nu))
    return alpha

def betaCallback(coe):
    # Beta is always 30 deg w/r to rsw frame
    return 30*np.pi/180

maneuver.addPerturbation("thrust")
maneuver.thrustProfile = (alphaCallback,betaCallback)
maneuver.propagate(60*60*24*10,60)
```
You can checkout many examples of the above in <a href="https://github.com/MrPapasFritas/frames-days/blob/master/Demo - Single Orbital Parameter Modification.ipynb">this demo</a>.

I've also implemented a method to define a target orbit and let the thrust orientation be calculated in order to reach that target orbit:
```python
# Add thrust
maneuvers.addPerturbation("thrust")
# Set Target Orbit
targetOrbitElements = [50000e3,0.01,10*np.pi/180,None,None]
maneuvers.setTargetOrbit(targetOrbitElements)
# Propagate for 18 days
maneuvers.propagate(60*60*24*18, 60)
```
An example of the above can be seen in <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Target Run.ipynb">this demo</a>  

In each propagation the data is saved in the histiry of the maneuver:
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
```
Some variables can also be calculated a posteriori, for example:
```python
maneuvers.calculateEclipseHours()
maneuvers.calculatePower()
```
Some examples of usage: 

<strong>Note: many of the examples below use the `ipyvolume` to create 3D plots. It is very recommended to try and install this library.</strong>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Deorbiting.ipynb">Demo - Deorbiting</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Inclination Change.ipynb">Demo - Inclination Change</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Perturbations.ipynb">Demo - Perturbations</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Eclipse.ipynb">Demo - Eclipse</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Energy.ipynb">Demo - Energy</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Relative Motion.ipynb">Demo - Relative Motion</a>
- <strong><a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Single Orbital Parameter Modification.ipynb">Demo - Single Orbital Parameter Modification</a></strong>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Demo - Target Run.ipynb">Demo - Target Run</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Validation/Validation - Perturbations.ipynb">Validation - Perturbations</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Validation/Validation - STK.ipynb">Validation - STK</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Demos/Calculations - Delta-V.ipynb">Calculations - Delta-V</a>

Some studies I attempted:
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Studies/Study - Station Keeping.ipynb">Study - Station Keeping</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Studies/Study - Deorbiting.ipynb">Study - Deorbiting</a>
- <a href="https://github.com/imricardoramos/orbital-mechanics/blob/master/Notebooks/Studies/Study - Orbiting bodies.ipynb">Study - Orbiting Bodies</a>
