"""
plot k value kKnoppWallFunction
plot omega value at wall boundary leeOmegaWallFunction
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
betaStar = 0.09  # constant of the model
beta1 = 0.075  # constant of the model
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


def omegaKnopp(ustar, y1):
    """
    Inputs:
    - ustar: friction velocity in m/s
    - y1: distance of the first cell center to nearest wall
    """
    d0 = getD0(ustar)
    knP = kn * ustar / nuF  # roughness Reynolds number
    omega = ustar / (np.sqrt(betaStar) * kappa * d0)
    limiter = 6 * nuF / (beta1 * y1**2)
    omega = np.where(omega < limiter, omega, limiter)
    return omega


foamTimes = os.popen("foamListTimes").read()
timeList = foamTimes.split("\n")[:-1]
timeArr = np.array([float(t) for t in timeList])
ntimes = len(timeList)

# friction velocity values
ustarArr = np.zeros(ntimes)
# k values at wall in simulation
kWall = np.zeros(ntimes)
# k values at first cell center in simulation
kFC = np.zeros(ntimes)
# k values at wall from Knopp formula
kKnopp = np.zeros(ntimes)

# omega values at wall in simulation
omWall = np.zeros(ntimes)
# omega values at first cell center in simulation
omFC = np.zeros(ntimes)
# omega values at wall from formula
omKnopp = np.zeros(ntimes)

Xmesh, Ymesh, Zmesh = rdf.readmesh("./")
y1 = Zmesh[0]  # first cell center distance to wall

for i, time in enumerate(timeList):
    # specific shear stress exerted by the wall
    tauW = rdf.readvector(
        "./", time, "wallShearStress", boundary="roughWall",
        verbose=False)[0, 0]
    ustarArr[i] = np.sqrt(np.abs(tauW))
    # turbulent kinetic energy
    kWall[i] = rdf.readscalar(
        "./", time, "k", boundary="roughWall",
        verbose=False)[0]
    kFC[i] = rdf.readscalar("./", time, "k", verbose=False)[0]
    # omega field
    omWall[i] = rdf.readscalar(
        "./", time, "omega", boundary="roughWall",
        verbose=False)[0]
    omFC[i] = rdf.readscalar("./", time, "omega", verbose=False)[0]

print(f"friction velocity: {ustarArr} m/s")

# k field at latest time step saved
kTend = rdf.readscalar("./", "latestTime", "k")
# omega field at latest time step saved
omTend = rdf.readscalar("./", "latestTime", "omega")

kKnopp[:] = kKnoppEisfeld(ustarArr)
omKnopp[:] = omegaKnopp(ustarArr, y1)

relErrOm = (omWall - omKnopp) / omKnopp
relErrK = (kWall - kKnopp) / kKnopp

fig = plt.figure(figsize=(14, 7))
gs = fig.add_gridspec(2, 4)

axK1 = fig.add_subplot(gs[:, 0])
axK1.plot(kTend, Zmesh, marker="x", color="#0072B2")
axK1.scatter(
    kWall[-1], 0, marker="o", color="firebrick", label="wall value")
axK1.scatter(
    kKnopp[-1], 0, marker="+", color="forestgreen", label="formula")
axK1.legend()
# zoom on near wall area
zmX1, zmX2, = axK1.get_xlim()
zmY1 = -0.002
zmY2 = 0.003
zmAxK = axK1.inset_axes(
    [0.7, 0.3, 0.47, 0.47],
    xlim=(zmX1, zmX2), ylim=(zmY1, zmY2))
zmAxK.plot(kTend, Zmesh, marker="x", color="#0072B2")
zmAxK.scatter(kWall[-1], 0, marker="o", color="firebrick")
zmAxK.scatter(kKnopp[-1], 0, marker="+", color="forestgreen")
zmAxK.grid()

axK2 = fig.add_subplot(gs[0, 1])
axK2.scatter(
    timeArr, kWall, marker="x",
    color="firebrick", label="wall value")
axK2.scatter(
    timeArr, kFC, marker="d",
    color="steelblue", label="first cell")
axK2.scatter(
    timeArr, kKnopp, marker="+",
    color="forestgreen", label="Knopp")

axErrK = fig.add_subplot(gs[1, 1])
axErrK.scatter(timeArr, relErrK, color="steelblue")

axOm1 = fig.add_subplot(gs[:, 2])
axOm1.plot(omTend, Zmesh, marker="x", color="#0072B2")
axOm1.scatter(omWall[-1], 0, color="firebrick")
axOm1.scatter(omKnopp[-1], 0, marker="+", color="forestgreen")
# zoom on near wall areazmX1 = omWall[-1] - 70
zmX1, zmX2, = axOm1.get_xlim()
zmY1 = -0.002
zmY2 = 0.003
zmAxOm1 = axOm1.inset_axes(
    [0.4, 0.3, 0.47, 0.47],
    xlim=(zmX1, zmX2), ylim=(zmY1, zmY2))
zmAxOm1.plot(omTend, Zmesh, marker="x", color="#0072B2")
zmAxOm1.scatter(omWall[-1], 0, color="firebrick")
zmAxOm1.scatter(omKnopp[-1], 0, color="forestgreen")
zmAxOm1.grid()

axOm2 = fig.add_subplot(gs[0, 3])
axOm2.scatter(
    timeArr, omWall, marker="x",
    color="firebrick", label="wall value")
axOm2.scatter(
    timeArr, omFC, marker="d",
    color="steelblue", label="first cell")
axOm2.scatter(
    timeArr, omKnopp, marker="+",
    color="forestgreen", label="Knopp")

axErrOm = fig.add_subplot(gs[1, 3])
axErrOm.scatter(timeArr, relErrOm, color="steelblue")

axK1.set_title(r"$k$ profile")
axK1.set_xlabel(r"$k\,[m^2s^{-2}]$")
axK1.set_ylabel(r"$z\,[m]$")
axK1.grid()

axK2.set_title(r"$k$(t) at wall BC")
axK2.set_ylabel(r"$k\,[s^{-1}]$")
axK2.set_xticklabels([])
axK2.grid()
axK2.legend()

axErrK.set_ylabel("relative error")
axErrK.set_xlabel("time s")
axErrK.grid()

axOm1.set_title(r"$\omega$ profile")
axOm1.set_xlabel(r"$\omega\,[s^{-1}]$")
axOm1.set_yticklabels([])
axOm1.grid()

axOm2.set_title(r"$\omega$(t) at wall BC")
axOm2.set_ylabel(r"$\omega\,[s^{-1}]$")
axOm2.set_xticklabels([])
axOm2.grid()
axOm2.legend()

axErrOm.set_ylabel("relative error")
axErrOm.set_xlabel("time s")
axErrOm.grid()

fig.tight_layout()

plt.show()
