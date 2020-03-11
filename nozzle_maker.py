#rocket-nozzle-tools
__author__ = "Ben Appleby"
__email__ = "ben.appleby@sky.com","b6040585@my.shu.ac.uk"
__copyright__ = "Copyright 2019, Ben Appleby"

"""Bell & Aerospike Nozzle Generator - Created by Benjamin Appleby @ Sheffield Hallam University [b6040585@my.shu.ac.uk]"""

import thermodynamics
import fluid
import contour
import geometry
import export
import plot
import boundary_condition_config as BCS


class Nozzle(fluid.ThermoFluids):
    """
    Nozzle object for producting boundary conditions and point geometry for a rocket nozzle of desired size and thrust. 

        Args:
            export_flow_domain: Boolean flag for exporting the flow domain to a file named nozzle.txt
    """
    def __init__(self,export_to_fluent = False): 
        target_thrust = input("Enter Desired Thrust (N): ")
        max_exit_radius = float(input("Enter Maximum (Exit) Radius (mm): "))
        try:
            length_factor = float(input("Enter Length Factor (1 = Full, 0.5 = Half etc: "))
        except:
            length_factor = 1
            
        if length_factor == 0:
            length_factor = 1

        self.set_variables()
        self.set_flow_conditions(float(target_thrust),max_exit_radius)
        self.exit_radius = max_exit_radius

        print("\n| Boundary Conditions: |")
        print("Chamber Pressure: "+str(self.p_total)+" [Pa]")
        print("Chamber Temperature: "+str(self.t_total)+" [K]")
        print("Obtained Thrust: "+str(self.thrust)+" [N]")
        print("Oxidizer Pressure: "+str(thermodynamics.calc_oxidizer_temp(self.p_total, self.t_total, BCS.oxidizer_temp)))

        self.coordinates = contour.plot_nozzle(self.a_ratio, max_exit_radius, self.mach, self.gamma, length_factor) #Area Ratio, Spike Radius, Mach Number and Gamma
        self.coordinates.sort(key = sortFirst, reverse = False)
        if export_to_fluent:
            export.export_flow_domain(self.coordinates,max_exit_radius,self.throat_radius) #Produces the flow domain geometry
        
    def plot_graph(self):
        """Plots a graph of the nozzle coordinates."""
        plot.plotGraph(self.coordinates)

def sortFirst(val):
    return val[0]


aerospikeNozzle = Nozzle(export_to_fluent=True)
