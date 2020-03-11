#rocket-nozzle-tools
__author__ = "Ben Appleby"
__email__ = "ben.appleby@sky.com","b6040585@my.shu.ac.uk"
__copyright__ = "Copyright 2019, Ben Appleby"

import math
import numpy as np

import thermodynamics


def make_rao_spline(x_curve2,y_curve2,attach_angle,x_exit,y_exit):
    theta_n = math.radians(attach_angle)
    raoMatrix1 = np.array([[y_curve2**2,y_curve2,1],[y_exit**2,y_exit,1],[2*y_curve2,1,0]])
    raoMatrix1 = np.linalg.inv(raoMatrix1)
    raoMatrix2 = np.array([[x_curve2],[x_exit],[1/math.tan(theta_n)]])

    return np.dot(raoMatrix1,raoMatrix2)

def calc_nozzle_length(R_throat,A_ratio,half_angle):
    nozzle_length = (R_throat*(math.sqrt(A_ratio)-1))/math.tan(math.radians(half_angle))
    return nozzle_length

def calc_chamber_length(throat_radius):
    return throat_radius/math.tan(math.radians(45))

def calc_throat_angle(throat_radius,exit_radius,nozzle_length,x_offset1,y_offset1): #Iterative Approach to calculating the nozzle-throat attachment angle
    xDiff = 100
    bestAttachAngle = 45

    for attach_angle_dec in range (2500,7500):
        attach_angle = attach_angle_dec/100
        theta_rads = math.radians(attach_angle)
        
        
        x = 0.382*throat_radius*math.sin(theta_rads)+x_offset1
        y = 0.382*throat_radius-abs(0.382*throat_radius*math.cos(theta_rads))+y_offset1
        x2 =0
        doIterate=True
        while x2<int(x) and doIterate==True:
            
            rao_coefficients = make_rao_spline(float(x),float(y),attach_angle,nozzle_length,exit_radius) #Gets Rao polynomial coefficients from the matrix
            rao_a = float(rao_coefficients[0])
            rao_b = float(rao_coefficients[1])
            rao_c = float(rao_coefficients[2])

            x2 = rao_a*(y**2)+rao_b*y-rao_c
            y+=0.01
            doIterate=False

        if(abs(x2-x)<xDiff):
            xDiff = abs(x2-x)
            bestAttachAngle  = attach_angle
            
    return bestAttachAngle

def calc_chamber_angle(throat_radius): # Iterative Approach to calculating the chamber-throat attachment angle
    dydxDiff = 100
    bestAttachAngle = 45

    for attach_angle_dec in range (1500,9000):
        attach_angle = attach_angle_dec/100
        theta_rads = math.radians(attach_angle)

        x_offset_chamber = 1.5*throat_radius*math.sin(theta_rads)
        y1 = 1.5*throat_radius-abs(1.5*throat_radius*math.cos(theta_rads))+throat_radius
        y2 = abs(1.5*throat_radius*math.cos(theta_rads))
        
        if abs(y2-y1)<dydxDiff:
            dydxDiff = abs(y2-y1)
            bestAttachAngle  = attach_angle
            
    return bestAttachAngle

def get_radius_mm(area):
    radius = math.sqrt(area/math.pi)*10**3
    return radius

def get_area_m(radius):
    area = math.pi*(radius*10**-3)**2
    return area

def calc_cowl_radius(exit_radius,throat_radius,spikeouterradius):
    cowlRadius = get_radius_mm(get_area_m(throat_radius)+get_area_m(spikeouterradius))
    return cowlRadius