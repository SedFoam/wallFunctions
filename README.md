wallFunctions
============

[![Release](https://img.shields.io/badge/release-0.1-blue.svg)](http://github.com/SedFoam/wallFunctions)
[![OpenFOAM v25xx](https://img.shields.io/badge/OpenFOAM-v25xx-brightgreen.svg)](https://openfoam.com/)
[![OpenFOAM v24xx](https://img.shields.io/badge/OpenFOAM-v24xx-brightgreen.svg)](https://openfoam.com/)

This repository provides **different wall functions** implementing **rough boundary conditions** for turbulent quantities **k (turbulent kinetic energy)** and **omega (specific dissipation rate)**. Developed at **LEGI**, this library extends OpenFOAM's native wall function capabilities to account for surface roughness effects in CFD simulations.

Available Models
----------------

This library includes the following rough wall boundary condition models:
   Model Name                     | Description                                                                |
 |--------------------------------|----------------------------------------------------------------------------|
 | `knoppkWallFunction`           | Knopp K wall function for rough walls.                                     |
 | `knoppOmegaWallFunction`       | Knopp Omega wall function for rough walls.                                 |
 | `fuhrmanOmegaWallFunction`     | Fuhrman Omega wall function for rough walls.                               |
 | `leeOmegaWallFunction`         | Lee Omega wall function with roughness corrections.                        |
 | `wilcoxOmegaWallFunction`      | Wilcox (2006) Omega wall function with roughness extensions.               |


Status
------

The roughWallFunctions library is in continuous development.

Pull requests are encouraged!

Installation
------------

```bash
cd $WM_PROJECT_USER_DIR
git clone https://github.com/sedfoam/wallFunctions
cd wallFunctions
./Allwmake
```

Usage
-----
Integrate into your OpenFOAM case:

Add the library to your controlDict:
```bash
libs ("libroughWallBC.so");
```

Use the desired boundary condition in your 0/ or constant/ boundary field files.

Developers
----------

  * Cyrille Bonamy
  * Julien Chauchat
  * Matthias Renaud
  * Antoine Mathieu
  * Eduard Puig Montella
  * Remi Chassagne
  * Maxime Kaczmarek
  * Tim Nagel
                              
Acknowledgements
----------------

OpenFOAM is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.

