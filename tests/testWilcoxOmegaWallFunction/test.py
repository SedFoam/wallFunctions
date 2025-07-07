"""
Test values taken by omega at wall boundary for wilcoxOmegaWallFunction
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- running test rough channel wilcoxOmegaWallFunction --- ")

success = True
tolOm = 1e-3

# equivalent roughness length
kn = 9.9e-3
nuF = 1e-6  # water kinematic viscosity
betaK = 0.09  # constant of the model
kappa = 0.41  # von Karmann constant


def getSr(knP, knLim=25.):
    """
    return Sr coefficient
    Inputs:
    - knP: ndarray, roughness Reynolds number
    - knLim: float, transition between hydraulic
    smooth and rough regimes
    """
    SrLow = (50 / knP)**2
    SrHigh = 100 / knP
    Sr = np.where(knP < knLim, SrLow, SrHigh)
    return Sr


def omegaWilcox(ustar):
    """
    Inputs:
    - ustar: friction velocity in m/s
    """
    knP = kn * ustar / nuF  # roughness Reynolds number
    Sr = getSr(knP)
    omega = Sr * ustar**2 / nuF
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
omWilcox = np.zeros(ntimes)

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


omWilcox[:] = omegaWilcox(ustarArr)

relErr = np.abs(omWall - omWilcox) / omWilcox

if np.any(relErr > tolOm):
    success = False
    print(
        "ERROR! omega values OpenFOAM and "
        + "computed from wallShearStress are not matching"
        + f"\nrelative error on omega: {relErr}")
else:
    print(f"omega value OK, maximum relative error {np.max(relErr)}")


assert success
