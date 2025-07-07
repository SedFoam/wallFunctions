"""

"""

import numpy as np
from abc import ABC, abstractmethod


def omegaVis(y1, nuF=1e-6, beta1=0.075):
    """
    viscous solution for omega
    Inputs:
    - y1: ndarray, distance from wall boundary
    - nuF: float, fluid kinematic viscosity
    - beta1: float, closure constant, 0.075
    """
    return 6 * nuF / (beta1 * y1**2)


def omegaLog(y1, k, betaStar=0.09, kappa=0.41):
    """
    viscous solution for omega
    Inputs:
    - y1: ndarray, distance from wall boundary
    - k: ndarray, specific turbulent kinetic energy
    - betaStar: float, closure constant, 0.09
    - kappa: float, von karman constant
    """
    return np.sqrt(k) / (np.sqrt(betaStar) * kappa * y1)


def omegaSmooth(
        y1,
        k,
        blending,
        nuF=1e-6,
        beta1=0.075,
        betaStar=0.09,
        kappa=0.41
):
    """
    Inputs:
    - y1: ndarray, distance from wall boundary
    - k: ndarray, specific turbulent kinetic energy
    - blending: str, blending function name
    - nuF: float, fluid kinematic viscosity
    - beta1: float, closure constant, 0.075
    - betaStar: float, closure constant, 0.09
    - kappa: float, von karman constant
    """
    F = blendingFunction.create(blending)
    omega = F.blendOmega(
        y1, k, nuF=1e-6, beta1=0.075, betaStar=0.09, kappa=0.41)
    return omega


def SrFuhrman(knP, knLim=5.):
    """
    return Sr coefficient from Fuhrman (2010)
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
    Sr = SrFuhrman(knP)
    omega = Sr * ustar**2 / nuF
    return omega


def getSr(knP, knLim=25.):
    """
    return Sr coefficient from Wilcox
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
    Sr = SrWilcox(knP)
    omega = Sr * ustar**2 / nuF
    return omega


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


def omegaLee(ustar, y1):
    """
    Inputs:
    - ustar: friction velocity in m/s
    - y1: distance of the first cell center to nearest wall
    """
    d0 = getD0(ustar)
    knP = kn * ustar / nuF  # roughness Reynolds number
    omega = ustar * np.log(1 + y1/d0) / (
        np.sqrt(betaStar) * kappa * d0)
    limiter = 6 * nuF / (betaStar * y1**2)
    omega = np.where(omega < limiter, omega, limiter)
    return omega


def kKnoppEisfeld(ustar):
    """
    boundary condition for k, from Knopp (2009)
    Inputs:
    - ustar: friction velocity in m/s
    """
    knP = kn * ustar / nuF  # roughness Reynolds number
    # specific turbulent kinetic energy
    k = ustar**2 / np.sqrt(betaStar)
    k *= np.where(knP < 90, knP / 90, 1)
    return k


# - - - - - BLENDING FUNCTIONS - - - - - #

class blendingFunction(ABC):

    blend_types = {}

    @classmethod
    def register_blend_type(cls, blend_type):
        def decorator(subclass):
            cls.blend_types[blend_type] = subclass
            return subclass
        return decorator

    @classmethod
    def create(cls, bFuncName):
        if bFuncName not in cls.blend_types:
            raise ValueError(
                "blending function not supported: " + bFuncName)
        return cls.blend_types[bFuncName]()


@blendingFunction.register_blend_type("max")
class max(blendingFunction):

    def blendOmega(
            self,
            y1,
            k,
            nuF=1e-6,
            beta1=0.075,
            betaStar=0.09,
            kappa=0.41
    ):
        omVis = omegaVis(y1, nuF=nuF, beta1=beta1)
        omLog = omegaLog(y1, k, betaStar=betaStar, kappa=kappa)
        return np.where(omVis > omLog, omVis, omLog)


@blendingFunction.register_blend_type("binomial2")
class binomial2(blendingFunction):

    def __init__(self):
        self.n = 2.

    def blendOmega(
            self,
            y1,
            k,
            nuF=1e-6,
            beta1=0.075,
            betaStar=0.09,
            kappa=0.41
    ):
        omVis = omegaVis(y1, nuF=nuF, beta1=beta1)
        omLog = omegaLog(y1, k, betaStar=betaStar, kappa=kappa)
        return (omVis**self.n + omLog**self.n)**(1/self.n)


@blendingFunction.register_blend_type("exponential")
class exponential(blendingFunction):

    pass


@blendingFunction.register_blend_type("tanh")
class tanh(blendingFunction):

    pass
