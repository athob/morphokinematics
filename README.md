# morphokinematics

Python2.7 routines for computing kinematics and morphological diagnostics of particles-made 3D systems. These routines go with the study [Thob et al. 2018](http://arxiv.org/abs/1811.01954) where such measurements have been computed for simulated galaxies. Please acknowledge by citing the latter if these routines have been significant for a research project that leads to a publication. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for usage.

### Prerequisites

You need to have a working Python2.7 installation with the following libraries already installed in their proposed version:

```
numpy 1.15.4
scipy 1.1.0
```

### Installing

To use the routines, simply download and copy the main script [morphokinematicsdiagnostics.py](morphokinematicsdiagnostics.py) to your working directory. Then execute it within your python environment via the following command:

```
execfile("morphokinematicsdiagnostics.py")
```

or

```
from morphokinematicsdiagnostics import *
```

### Functions

This script contains 4 functions, and each are associated with a detailed help page that can be generated via the Python built-in ```help``` function.

#### cumsummedian

Computes and returns the weighted median value of an array ```a``` of some quantity, given an array ```weights``` of corresponding weights. It uses the following syntax:

```
cumsummedian(a,weights)
```

```a``` and ```weights``` must both be 1-dimensional arrays of same length.

#### center_and_unloop

Computes and returns the centered and unlooped coordinates of a collection of particles positions ```XYZ```, given an array ```XYZ0``` of the center coordinates and a quantity ```BoxL``` for the size of the looping box. It uses the following syntax:

```
center_and_unloop(XYZ,XYZ0,BoxL)
```

```XYZ``` must be a 2-dimensional (n, 3) array with the 3 coordinates along the 2nd dimension, while ```XYZ0``` must contain the centre 3 coordinates as an array_like of length 3.

#### kinematics_diagnostics

Computes and returns the kinematics diagnostics and properties of a collection of particles, given their positions ```XYZ```, their masses ```mass```, their velocities ```Vxyz``` and their specific binding energies ```PBE```, within an aperture radius ```aperture```. It uses the following syntax:

```
kinematics_diagnostics(XYZ,mass,Vxyz,PBE,aperture)
```

```XYZ``` and ```Vxyz``` must be 2-dimensional (n, 3) arrays with the 3 coordinates along the 2nd dimension, of same lengths with ```mass``` and ```PBE``` that must be 1-dimensional arrays.

#### morphological_diagnostics

Computes and returns the morphological diagnostics and properties of a collection of particles, given their positions ```XYZ```, their masses ```mass``` and their velocities ```Vxyz```, within an aperture radius ```aperture```. It uses the following syntax:

```
morphological_diagnostics(XYZ,mass,Vxyz,aperture=0.03)
```

```XYZ``` and ```Vxyz``` must be 2-dimensional (n, 3) arrays with the 3 coordinates along the 2nd dimension, of same lengths with ```mass``` that must be a 1-dimensional array.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## Authors

* **Adrien C.R. Thob** - *Initial work* - [athob](http://github.com/athob)

See also the list of [contributors](http://github.com/athob/morphokinematics/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The owner gratefully acknowledges the support of the [Royal Society](http://royalsociety.org/) via a doctoral studentship and the [Liverpool John Moores University](http://www.ljmu.ac.uk/) - [Astrophysics Research Institute](http://www.astro.ljmu.ac.uk) where this doctorate degree has been carried out
* These routines have been contributed by insightful guidance of the owner doctoral supervisory team, [Robert A. Crain](http://www.astro.ljmu.ac.uk/~astrcrai/), [Ian G. McCarthy](http://www.astro.ljmu.ac.uk/~igm/) and [Philip A. James](http://www.ljmu.ac.uk/about-us/staff-profiles/faculty-of-engineering-and-technology/astrophysics-research-institute/philip-james)
* Thanks to all the other co-authors of the study [Thob et al. 2018](http://arxiv.org/abs/1811.01954) for which these routines were initially developed
* Thanks to [Simon Pfeifer](https://github.com/SimonPfeifer) for his thorough proof testing and his improvements suggestions
