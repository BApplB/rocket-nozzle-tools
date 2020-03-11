#rocket-nozzle-tools
__author__ = "Ben Appleby"
__email__ = "ben.appleby@sky.com","b6040585@my.shu.ac.uk"
__copyright__ = "Copyright 2019, Ben Appleby"

import math

import thermodynamics
import geometry

def plot_nozzle(area_ratio,exit_radius,mach,gamma,lengthfactor,is_aerospike = True):
    exit_angle = math.degrees(thermodynamics.calc_prandtl_meyer_angle(mach,gamma)*mach) #Gets the Prandtl-Meyer angle to reduce losses

    area_ratio = area_ratio
    calc_exit_area = geometry.get_area_m(exit_radius)
    throat_area = calc_exit_area/area_ratio
    throat_radius = geometry.get_radius_mm(throat_area)
    chamber_length = geometry.calc_chamber_length(throat_radius)
    nozzle_length = geometry.calc_nozzle_length(throat_radius,area_ratio,exit_angle)
    truncation_position = lengthfactor*nozzle_length

    if not is_aerospike:
        nozzle_length = truncation_position
    

    x_offset1 = 1.5*throat_radius*math.sin(0)
    y_offset1 = throat_radius

    chamber_attach_angle = geometry.calc_chamber_angle(throat_radius)              
    throat_attach_angle = geometry.calc_throat_angle(throat_radius,exit_radius,nozzle_length,x_offset1,y_offset1)
    

    x_offset_chamber = 1.5*throat_radius*math.sin(-math.radians(chamber_attach_angle))
    y_radius_chamber = 1.5*throat_radius
    
    #Generating the points comprising the aerospike

    combustion_chamber =  make_combustion_chamber(chamber_length,throat_radius,y_radius_chamber,x_offset_chamber)
    combustion_chamber_curve =  make_combustion_chamber_curve(chamber_attach_angle,throat_radius,x_offset_chamber)
    converging_throat_curve = make_converging_throat_curve(chamber_attach_angle,throat_radius)
    diverging_throat_curve = make_diverging_throat_curve(throat_attach_angle,throat_radius,x_offset1,y_offset1)

    coordinates = combustion_chamber + combustion_chamber_curve + converging_throat_curve + diverging_throat_curve
    
    rao_coefficients = geometry.make_rao_spline(float(coordinates[-1][0]),float(coordinates[-1][1]),throat_attach_angle,nozzle_length,exit_radius)
    rao_a = float(rao_coefficients[0])
    rao_b = float(rao_coefficients[1])
    rao_c = float(rao_coefficients[2])

    x, y = coordinates[-1] #Set x and y to the last coordinates added

    while (x>=int(coordinates[-1][0]) and x<=int(truncation_position) and y<exit_radius): #Rao Thrust Curve
        x = rao_a*(y**2)+rao_b*y+rao_c
        coordinates.append([x,y])
        y += 0.05

    if is_aerospike: #If the create aerospike toggle is enabled
        coordinates = coordinates_as_aerospike(coordinates,exit_radius,truncation_position,nozzle_length)
        

    print("\nAREA RATIO: "+str(area_ratio)+" \nMACH: "+str(mach))
    print("A = "+str(calc_exit_area))
    print("A* = "+str(throat_area))
    print("R* = "+str(throat_radius))
    print("Ln = "+str(nozzle_length))
    return coordinates



def make_combustion_chamber(chamber_length,throat_radius,y_radius_chamber,x_offset_chamber):
    coordinates = []
    for L in range(0,int(chamber_length*10)):
        x = (2*x_offset_chamber-1.5*throat_radius*math.sin(0)) - (L/10)
        y = y_radius_chamber
        coordinates.append([x,y])
    return coordinates

def make_combustion_chamber_curve(chamber_attach_angle,throat_radius,x_offset_chamber):
    coordinates = []
    for theta in range(-int(chamber_attach_angle*10),0): #Combustion Chamber Curve
        thetarads = 0.1*theta*(math.pi/180)
        x = 2*x_offset_chamber-1.5*throat_radius*math.sin(thetarads)
        y = abs(1.5*throat_radius*math.cos(thetarads))
        coordinates.append([x,y])
    return coordinates

def make_converging_throat_curve(chamber_attach_angle,throat_radius):
    coordinates = []
    for theta in range(-int(chamber_attach_angle*10),0): #Converging Curve
        thetarads = 0.1*theta*(math.pi/180)
        x = 1.5*throat_radius*math.sin(thetarads)
        y = 1.5*throat_radius-abs(1.5*throat_radius*math.cos(thetarads))+throat_radius
        coordinates.append([x,y])
    return coordinates

def make_diverging_throat_curve(throat_attach_angle,throat_radius,x_offset1,y_offset1):
    coordinates = []
    for theta in range(0,int(throat_attach_angle*10)): #Diverging Curve
        thetarads = 0.1*theta*(math.pi/180)
        x = 0.382*throat_radius*math.sin(thetarads)+x_offset1
        y = 0.382*throat_radius-abs(0.382*throat_radius*math.cos(thetarads))+y_offset1
        coordinates.append([x,y])
    return coordinates

def coordinates_as_aerospike(coordinates,exit_radius,truncation_position,nozzle_length):
    index = 0
    while index < len(coordinates):
        coordinates[index][1] = abs(coordinates[index][1]-exit_radius) #Bring the nozzle wall exit to y=0
        index += 1

    if(truncation_position!=nozzle_length):
        coordinates.append([truncation_position,0])
    return coordinates
