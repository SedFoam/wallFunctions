"""
Test values taken by omega at wall boundary for leeOmegaWallFunction
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- running test rough channel leeOmegaWallFunction --- ")

success = True
tolOm = 1e-4

# equivalent roughness length
kn = 9.9e-3
nuF = 1e-6  # water kinematic viscosity
betaStar = 0.09  # constant of the model
beta1 = 0.075  # turbulence model closure constant
kappa = 0.41  # von Karmann constant


def getD0(ustar):
    """
    corrected hydrodynamical roughness, Cheng-Hsien Lee (2017)
    Inputs:
    - ustar: ndarray, friction velocity"""
    knP = kn * ustar / nuF  # roughness Reynolds number
    d0 = 0.03 * kn
    d0 *= np.where(knP < 30, (knP/30)**(2/3), 1)
    d0 *= np.where(knP < 45, (knP/45)**(1/4), 1)
    d0 *= np.where(knP < 60, (knP/60)**(1/4), 1)
    return d0


def kKnoppEisfeld(ustar):
    """
    Inputs:
    - ustar: friction velocity in m/s
    """
    knP = kn * ustar / nuF  # roughness Reynolds number
    # specific turbulent kinetic energy
    k = np.where(knP < 90, knP / 90, 1) * ustar**2 / np.sqrt(betaStar)
    return k


def omegaLee(ustar, y1):
    """
    Inputs:
    - ustar: friction velocity in m/s
    - y1: distance of the first cell center to nearest wall
    """
    d0 = getD0(ustar)
    knP = kn * ustar / nuF  # roughness Reynolds number
    omega = ustar * np.log(1 + y1/d0) / (np.sqrt(betaStar) * kappa * d0)
    limiter = 6 * nuF / (beta1 * y1**2)
    omega = np.where(omega < limiter, omega, limiter)
    return omega


foamTimes = os.popen("foamListTimes").read()
timeList = foamTimes.split("\n")[:-1]
timeArr = np.array([float(t) for t in timeList])
ntimes = len(timeList)

# friction velocity values
ustarArr = np.zeros(ntimes)
# omega values at wall in simulation
omWall = np.zeros(ntimes)
# omega values at first cell center in simulation
omFC = np.zeros(ntimes)
# omega values at wall from formula
omLee = np.zeros(ntimes)

Xmesh, Ymesh, Zmesh = rdf.readmesh("./", verbose=False)
y1 = Zmesh[0]  # first cell center distance to wall

for i, time in enumerate(timeList):
    # specific shear stress exerted by the wall
    tauW = rdf.readvector(
        "./", time, "wallShearStress",
        boundary="roughWall", verbose=False)[0, 0]
    ustarArr[i] = np.sqrt(np.abs(tauW))
    omWall[i] = rdf.readscalar(
        "./", time, "omega", boundary="roughWall",
        verbose=False)[0]
    omFC[i] = rdf.readscalar(
        "./", time, "omega",
        verbose=False)[0]


omLee[:] = omegaLee(ustarArr, y1)

relErr = np.abs(omWall - omLee) / omLee

if np.any(relErr > tolOm):
    success = False
    print(
        "ERROR! omega values OpenFOAM and "
        + "computed from wallShearStress are not matching"
        + f"\nrelative error on omega: {relErr}")
else:
    print(f"omega value OK, maximum relative error {np.max(relErr)}")


assert success
