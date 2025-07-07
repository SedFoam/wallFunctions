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


def nutLee(ustar, y1):
    """
    turbulent eddy viscosity at wall boundary from Cheng-Hsien Lee
    Inputs:
    - ustar: friction velocity in m/s
    - y1: distance of the first cell center to nearest wall
    """
    d0 = getD0(ustar)
    nut = kappa * ustar * y1 / np.log(1 + y1/d0)
    return nut


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
omLee = np.zeros(ntimes)

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
# turbulent eddy viscosity field
nutField = rdf.readscalar("./", "latestTime", "nut")
nutWall = rdf.readscalar(
    "./", "latestTime", "nut", boundary="roughWall", verbose=False)[0]

kKnopp[:] = kKnoppEisfeld(ustarArr)
omLee[:] = omegaLee(ustarArr, y1)

relErr = (omWall - omLee) / omLee

fig = plt.figure(figsize=(15, 7))
gs = fig.add_gridspec(2, 4)

axK = fig.add_subplot(gs[:, 0])
axK.plot(
    kTend, Zmesh, marker="x", color="steelblue")
axK.scatter(
    kWall[-1], 0, marker="o", color="firebrick", label="wall value")
axK.scatter(
    kKnopp[-1], 0, marker="+", color="forestgreen", label="formula")
axK.legend()

# zoom on near wall area
zmX1, zmX2, = axK.get_xlim()
zmY1 = -0.002
zmY2 = 0.003
zmAxK = axK.inset_axes(
    [0.7, 0.3, 0.47, 0.47],
    xlim=(zmX1, zmX2), ylim=(zmY1, zmY2))
zmAxK.plot(kTend, Zmesh, marker="x", color="steelblue")
zmAxK.scatter(kWall[-1], 0, marker="o", color="firebrick")
zmAxK.scatter(kKnopp[-1], 0, marker="+", color="forestgreen")
zmAxK.grid()

axOm1 = fig.add_subplot(gs[:, 1])
axOm1.plot(omTend, Zmesh, marker="x", color="#0072B2")
axOm1.scatter(omWall[-1], 0, marker="o", color="firebrick")
axOm1.scatter(omLee[-1], 0, marker="+", color="forestgreen")

zmX1, zmX2, = axOm1.get_xlim()
zmY1 = -0.002
zmY2 = 0.003
zmAxOm1 = axOm1.inset_axes(
    [0.4, 0.3, 0.47, 0.47],
    xlim=(zmX1, zmX2), ylim=(zmY1, zmY2))
zmAxOm1.plot(omTend, Zmesh, marker="x", color="steelblue")
zmAxOm1.scatter(omWall[-1], 0, marker="o", color="firebrick")
zmAxOm1.scatter(omLee[-1], 0, marker="+", color="forestgreen")
zmAxOm1.grid()

axOm2 = fig.add_subplot(gs[0, 2])
axOm2.scatter(
    timeArr, omWall, marker="o",
    color="firebrick", label="wall value")
axOm2.scatter(
    timeArr, omFC, marker="x",
    color="steelblue", label="first cell")
axOm2.scatter(
    timeArr, omLee, marker="+",
    color="forestgreen", label="Cheng-Hsien Lee")

axErr = fig.add_subplot(gs[1, 2])
axErr.scatter(timeArr, relErr, color="steelblue")

axNut = fig.add_subplot(gs[:, 3])
axNut.plot(nutField, Zmesh, marker="x", color="steelblue")
axNut.scatter(nutWall, 0, marker="o", color="firebrick")
axNut.scatter(nutLee(ustarArr[-1], y1), 0,  marker="+", color="forestgreen")
axNut.axline((0, 0), slope=1/(ustarArr[-1]*kappa), color="black", ls="dashed")

axK.set_title(r"$k$ profile")
axK.set_xlabel(r"$k\,[m^2s^{-2}]$")
axK.set_ylabel(r"$z\,[m]$")
axK.grid()

axOm1.set_title(r"$\omega$ profile")
axOm1.set_xlabel(r"$\omega\,[s^{-1}]$")
axOm1.set_yticklabels([])
axOm1.grid()

axOm2.set_title(r"$\omega$(t) at wall BC")
axOm2.set_ylabel(r"$\omega\,[s^{-1}]$")
axOm2.set_xticklabels([])
axOm2.grid()
axOm2.legend()

axErr.set_ylabel("relative error")
axErr.set_xlabel("time s")
axErr.grid()

axNut.set_xlabel(r"$\nu_t\,[m^2.s^{-1}]$")
axNut.set_ylabel(r"$z\,[m]$")
axNut.grid()

fig.tight_layout()

plt.show()
