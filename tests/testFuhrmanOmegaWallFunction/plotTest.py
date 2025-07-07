"""
plot omega value at wall boundary fuhrmanOmegaWallFunction
compare result from OpenFOAM and direct computation from shear stress
"""


import numpy as np
from fluidfoam import readof as rdf
import matplotlib.pyplot as plt
import os


plt.rcParams["font.size"] = 12

# equivalent roughness length
kn = 9.9e-3
nuF = 1e-6  # water kinematic viscosity
betaK = 0.09  # constant of the model
kappa = 0.41  # von Karmann constant


def getSr(knP, knLim=5.):
    """
    return Sr coefficient
    Inputs:
    - knP: ndarray, roughness Reynolds number
    - knLim: float, transition between hydraulic
    smooth and rough regimes
    """
    SrLow = (200 / knP)**2
    SrHigh = (100 / knP) + (
        (200 / knP)**2 - (100 / knP)) * np.exp(knLim-knP)
    Sr = np.where(knP < knLim, SrLow, SrHigh)
    return Sr


def omegaFuhrman(ustar):
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
omFuhrman = np.zeros(ntimes)

Xmesh, Ymesh, Zmesh = rdf.readmesh("./")
y1 = Zmesh[0]  # first cell center distance to wall

for i, time in enumerate(timeList):
    # specific shear stress exerted by the wall
    tauW = rdf.readvector(
        "./", time, "wallShearStress",
        boundary="roughWall", verbose=False)[0, 0]
    ustarArr[i] = np.sqrt(np.abs(tauW))
    omWall[i] = rdf.readscalar(
        "./", time, "omega", boundary="roughWall", verbose=False)[0]
    omFC[i] = rdf.readscalar("./", time, "omega")[0]

print(f"friction velocity: {ustarArr} m/s")

# omega field at latest time step saved
omTend = rdf.readscalar("./", "latestTime", "omega")

omFuhrman[:] = omegaFuhrman(ustarArr)

relErr = np.abs(omWall - omFuhrman) / omFuhrman


fig = plt.figure(figsize=(10, 7))
gs = fig.add_gridspec(2, 2)

axOm1 = fig.add_subplot(gs[:, 0])
axOm1.plot(omTend, Zmesh, marker="x", color="#0072B2")
axOm1.scatter(omWall[-1], 0, color="firebrick")
# zoom on near wall areazmX1 = omWall[-1] - 70
zmX1 = 100
zmX2 = omWall[-1] + 20
zmY1 = -0.002
zmY2 = 0.003
zmAx = axOm1.inset_axes(
    [0.4, 0.5, 0.47, 0.47],
    xlim=(zmX1, zmX2), ylim=(zmY1, zmY2))
zmAx.plot(omTend, Zmesh, marker="x", color="#0072B2")
zmAx.scatter(omWall[-1], 0, color="firebrick")
zmAx.grid()

axOm2 = fig.add_subplot(gs[0, 1])
axOm2.scatter(
    timeArr, omWall, marker="o",
    color="firebrick", label="wall value")
axOm2.scatter(
    timeArr, omFC, marker="x",
    color="steelblue", label="first cell")
axOm2.scatter(
    timeArr, omFuhrman, marker="+",
    color="forestgreen", label="Fuhrman")

axErr = fig.add_subplot(gs[1, 1])
axErr.scatter(timeArr, relErr, color="steelblue")

axOm1.set_title(r"$\omega$ profile")
axOm1.set_xlabel(r"$\omega\,[s^{-1}]$")
axOm1.set_ylabel(r"$z\,[m]$")
axOm1.grid()

axOm2.set_title(r"$\omega$(t) at wall BC")
axOm2.set_ylabel(r"$\omega\,[s^{-1}]$")
axOm2.set_xticklabels([])
axOm2.grid()
axOm2.legend()

axErr.set_ylabel("relative error")
axErr.set_xlabel("time s")
axErr.grid()

fig.tight_layout()

plt.show()
