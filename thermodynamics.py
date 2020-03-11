#rocket-nozzle-tools
__author__ = "Ben Appleby"
__email__ = "ben.appleby@sky.com","b6040585@my.shu.ac.uk"
__copyright__ = "Copyright 2019, Ben Appleby"

import math

"""Common thermodynamics functions"""

def x_function(mach,gamma):
    gp1 = gamma+1
    gm1 = gamma-1
    return mach*math.sqrt(gamma)*(1+(gm1/2)*mach**2)**(-gp1/(2*gm1))

def pratio(mach,gamma):
    pressure_ratio = (1+((gamma-1)/2)*mach**2)**(gamma/(gamma-1))
    return pressure_ratio

def aratio(mach,gamma):
    area_ratio = ((gamma+1)/2)**(-((gamma+1)/(gamma-1)/2)) / mach * (1 + mach**2 * (gamma-1)/2)**((gamma+1)/(gamma-1)/2)
    return area_ratio

def rhoratio(mach,gamma):
    rho_ratio = (1+((gamma-1)/2)*mach**2)**(-1/(gamma-1))
    return rho_ratio

def tratio(mach,gamma):
    temperature_ratio = (1+((gamma-1)/2)*mach**2)**(-1)
    return temperature_ratio

def pcrit(gamma):
    pressure_critical = (2/(gamma+1))**(gamma/(gamma-1))
    return pressure_critical

def calc_oxidizer_temp(chamberPressure,chamberTemp,OxTemp):
    return (chamberPressure*OxTemp)/chamberTemp

def calc_prandtl_meyer_angle(mach,gamma):
    calc_pm = math.sqrt(gamma+1/gamma-1)*math.atan(math.sqrt((gamma-1/gamma+1)*(mach**2-1)))-math.atan(math.sqrt(mach**2-1))
    return calc_pm

def sound_speed(gamma,specificgasconstant,temperature):
    calc_sos = math.sqrt(gamma*specificgasconstant*temperature)
    return calc_sos