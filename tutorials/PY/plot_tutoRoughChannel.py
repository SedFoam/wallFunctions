"""
smooth Channel case from Fuhrman (2010)
compare OpenFoam results with experimental measurements
use of omegaWallFunction on bottom wall boundary
"""

import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readof as rdf
import os

plt.rcParams["font.size"] = 12


pathData = "../DATA/Fuhrman2010/"
time = "latestTime"


caseList = [
    {"path": "../roughChannel/Wilcox",
     "time": "latestTime",
     "label": "Wilcox",
     "ls": "solid",
     "color": "#0072B2"},
    {"path": "../roughChannel/Fuhrman",
     "time": "latestTime",
     "label": "Fuhrman",
     "ls": "dashed",
     "color": "#D55E00"},
    {"path": "../roughChannel/Knopp",
     "time": "latestTime",
     "label": "Knopp",
     "ls": "dashdot",
     "color": "#009E73"},
    {"path": "../roughChannel/Lee",
     "time": "latestTime",
     "label": "Lee",
     "ls": "dotted",
     "color": "#CC79A7"}]

# physical and experimental parameters/data
ustarExp = 0.021  # experimental friction velocity
nuF = 9.6e-7  # water kinematic viscosity
Hwater = 6.7e-2  # water height
ksExp = 9.9e-3  # Nikuradse equivalent roughness height

# load experimental data
# u/ustar, z/ks
uxDat, zuksDat = np.loadtxt(
    pathData + "u_rough.txt", unpack=True, delimiter=" ")
# k/ustar^2, z/ks
kDat, zkksDat = np.loadtxt(
    pathData + "k_rough.txt", unpack=True, delimiter=" ")


fig1, (axU, axK, axOm, axNut) = plt.subplots(
    ncols=4, figsize=(12, 6))

# zoom on omega ax
zmX1 = 50
zmX2 = 500
zmY1 = -0.002
zmY2 = 0.05
zmAx = axOm.inset_axes(
    [0.6, 0.5, 0.47, 0.47],
    xlim=(zmX1, zmX2), ylim=(zmY1, zmY2))

fig2, (axUstar, axYplus, axKsPlus) = plt.subplots(
    nrows=3, figsize=(8, 8))

for case in caseList:
    pathCase = case["path"]
    color = case["color"]
    ls = case["ls"]
    label = case["label"]

    lw = 2.
    if ls == "dotted":
        lw = 3.

    Xmesh, Ymesh, Zmesh = rdf.readmesh(pathCase)

    # velocity field, x-direction
    Ux = rdf.readvector(pathCase, time, "U")[0]
    axU.plot(
        Ux, Zmesh/Hwater, lw=lw, ls=ls, color=color, label=label)

    # k, specific turbulent kinetic enery
    kField = rdf.readscalar(pathCase, time, "k")
    axK.plot(
        kField, Zmesh/Hwater, lw=lw, ls=ls, color=color)
    # omega
    omField = rdf.readscalar(pathCase, time, "omega")
    axOm.plot(
        omField, Zmesh/Hwater, lw=lw, ls=ls, color=color)
    zmAx.plot(
        omField, Zmesh/Hwater, lw=lw, ls=ls, color=color)

    # turbulent eddy viscosity
    nutField = rdf.readscalar(pathCase, time, "nut")
    axNut.plot(
        nutField, Zmesh/Hwater, lw=lw, ls=ls, color=color)

    foamTimes = os.popen(f"foamListTimes -case {pathCase}").read()
    timeList = foamTimes.split("\n")[:-1]
    timeArr = np.array([float(t) for t in timeList])
    ntimes = len(timeList)
    ustarArr = np.zeros(ntimes)
    yPlusArr = np.zeros(ntimes)

    for i, t in enumerate(timeList):
        # read wall shear stress exerted by wall
        tauWall = rdf.readvector(
            pathCase, t, "wallShearStress",
            boundary="roughWall", verbose=False)[0, 0]
        ustarArr[i] = np.sqrt(np.abs(tauWall))
        yPlusArr[i] = rdf.readscalar(
            pathCase, t, "yPlus", boundary="roughWall", verbose=False)[0]
    ksPlusArr = ksExp * ustarArr / nuF
    axUstar.plot(
        timeArr, ustarArr, lw=lw, ls=ls, color=color, label=label)
    axYplus.plot(
        timeArr, yPlusArr, lw=lw, ls=ls, color=color)
    axKsPlus.plot(
        timeArr, ksPlusArr, lw=lw, ls=ls, color=color)

axU.scatter(
    uxDat*ustarExp, zuksDat*ksExp / Hwater,
    color="firebrick", marker="x", label="Fuhrman (2010)")
axU.set_xlabel(r"$u_x\,[m.s^{-1}]$")
axU.set_ylabel(r"$z/H$")
axUstar.legend(loc="upper left")
axU.grid()
axU.legend()

axK.scatter(
    kDat*ustarExp**2, zkksDat*ksExp / Hwater,
    color="firebrick", marker="x")
axK.set_xlabel(r"$k\,[m^2.s^{-2}]$")
axK.set_yticklabels([])
axK.grid()

axOm.set_xlabel(r"$\omega$")
axOm.set_xscale("log")
axOm.set_yticklabels([])
axOm.grid()
zmAx.grid()

axNut.set_xlabel(r"$\nu_t$")
axNut.set_yticklabels([])
axNut.grid()

axUstar.set_ylabel(r"$u_*\,[m.s^{-1}]$")
axUstar.axhline(ustarExp, color="black", ls="dashed", label=r"$u_*$ exp")
axUstar.set_xticklabels([])
axUstar.grid()

axYplus.set_ylabel(r"$z^+$")
axYplus.set_xticklabels([])
axYplus.grid()

axKsPlus.set_xlabel("t [s]")
axKsPlus.set_ylabel(r"$k_s^+$")
axKsPlus.grid()


fig1.tight_layout()

fig2.tight_layout()

plt.show()
