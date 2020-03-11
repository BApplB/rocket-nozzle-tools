#rocket-nozzle-tools
__author__ = "Ben Appleby"
__email__ = "ben.appleby@sky.com","b6040585@my.shu.ac.uk"
__copyright__ = "Copyright 2019, Ben Appleby"

import math

import thermodynamics
import geometry
import boundary_condition_config as BCS

"""Thermofluid Dynamics module for determining nozzle flow conditions."""

class ThermoFluids:
    def __init__(self):
        self.set_variables()
    
    def set_variables(
        self,
        inlet_specifictemperature = BCS.inlet_specifictemperature, outlet_staticpressure = BCS.outlet_staticpressure,
        gamma = BCS.gamma, specific_gas_constant = BCS.specific_gas_constant,
        inlet_adiabatictemp = BCS.inlet_adiabatictemp ,oxidiser_fraction = BCS.oxidiser_fraction
    ):
        self.inlet_specifictemperature = inlet_specifictemperature
        self.outlet_staticpressure = outlet_staticpressure
        self.gamma = gamma
        self.specific_gas_constant = specific_gas_constant
        self.inlet_adiabatictemp = inlet_adiabatictemp
        self.oxidiser_fraction = oxidiser_fraction
    
    def set_flow_conditions(self,thrust,max_exit_radius):
        self.p_total,self.mach,self.thrust,self.a_ratio,self.t_total,self.throat_radius = self.calc_flow_conditions(thrust,max_exit_radius)

    def calc_flow_conditions(self,thrust,max_exit_radius): # Function for determining the chamber boundary conditions
        gp1 = self.gamma+1
        gm1 = self.gamma-1
        temp = self.inlet_specifictemperature
        mach = 4
        calc_exit_area = geometry.get_area_m(max_exit_radius)
        for iteration in range(0,10000):

            pressure_ratio = thermodynamics.pratio(mach,self.gamma)
            pressure_critical = thermodynamics.pcrit(self.gamma)
            area_ratio = thermodynamics.aratio(mach,self.gamma)
            rho_ratio = thermodynamics.rhoratio(mach,self.gamma)
            temperature_ratio = thermodynamics.tratio(mach,self.gamma)

            exit_temperature = temperature_ratio*temp
            
            exit_velocity = mach*thermodynamics.sound_speed(self.gamma,self.specific_gas_constant,exit_temperature)
            throat_velocity = mach*thermodynamics.sound_speed(self.gamma,self.specific_gas_constant,temp)

            pressure_total = (self.outlet_staticpressure)*pressure_ratio
            throat_pressure = pressure_total/pressure_critical
            chamber_density = pressure_total/(self.specific_gas_constant*temp)
            
            throat_area = calc_exit_area/area_ratio
            
            throat_radius = geometry.get_radius_mm(throat_area)
            
            massflowrate = ((pressure_total*throat_area)/math.sqrt(self.specific_gas_constant*temp))*thermodynamics.x_function(mach,self.gamma)

            calc_oxidiser_massflowrate = massflowrate*self.oxidiser_fraction
            temp = self.inlet_specifictemperature*calc_oxidiser_massflowrate #Calculates the adiabatic flame temperature based on oxidiser mass flow

            if (temp>self.inlet_adiabatictemp): #Prevents the temperature from exceeding the limits of the model
                temp = self.inlet_adiabatictemp
                if (iteration==9999):
                    print("\n**WARNING: Temperature Exceeded Adiabatic Flame Temp!**\n")
            
            calc_thrust = massflowrate*exit_velocity
            
            if calc_thrust > thrust:
                mach +=0.01
            elif calc_thrust < thrust:
                mach -=0.01
            else:
                return [pressure_total,mach]
        print("\n--Debug--")
        print("m_dot: "+str(massflowrate)+" [kg/s]")
        print("rho0: "+str(chamber_density)+" [kg/m^3]")
        print("P*_crit: "+str(throat_pressure)+" [Pa]")
        print("P_ratio: "+str(pressure_ratio))
        print("R_throat: "+str(throat_radius))
        print("---------")
        return [pressure_total,mach,calc_thrust,area_ratio,temp,throat_radius]