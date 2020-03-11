#rocket-nozzle-tools
__author__ = "Ben Appleby"
__email__ = "ben.appleby@sky.com","b6040585@my.shu.ac.uk"
__copyright__ = "Copyright 2019, Ben Appleby"

"""Functions for exporting the geometry as a flow domain."""

import geometry

FLOW_DOMAIN_MULTIPLIER = 5

def export_flow_domain(coordinates,max_exit_radius,throat_radius,is_aerospike=False):
    file = open("nozzle.txt",'w')

    chambery = coordinates[0][1]
    chamberx = coordinates[0][0]

    maxy = max_exit_radius
    max_spikeradius = 0
    maxx = 0
    xstep = coordinates[int(len(coordinates)-1)][0] - coordinates[int(len(coordinates)-2)][0]
    point = 0
    pointer = 1

    while point in range(0,len(coordinates)):
        file.write("1"+"    "+str(pointer)+"   "+str(coordinates[point][0])+"    "+str(coordinates[point][1])+"   "+"0\n")
        if float(coordinates[point][1])>max_spikeradius:
           max_spikeradius = float(coordinates[point][1])
        if float(coordinates[point][1])>maxy:
           maxy = float(coordinates[point][1])
        if float(coordinates[point][0])>maxx:
           maxx = float(coordinates[point][0])
        point +=2
        pointer+=1
        
    maxx += 0.00001
        
    #Used to produce the rest of the flow domain
    if is_aerospike:
        maxy = geometry.calc_cowl_radius(max_exit_radius,throat_radius,max_spikeradius)
        file.write("1"+"    "+str(pointer)+"   "+str(maxx+xstep)+"    "+"0"+"   "+"0\n") #Terminates the loop
        
        domainMatrix = [(2,1,maxx+xstep,0,0),(2,2,chamberx,0,0),(3,1,chamberx,0,0),(3,2,chamberx,chambery,0),(4,1,chamberx,chambery,0),(4,2,chamberx,maxy,0),
                        (5,1,chamberx,maxy,0),(5,2,0,maxy,0),(6,1,0,maxy,0),(6,2,0,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),(7,1,0,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),(7,2,maxx*12*FLOW_DOMAIN_MULTIPLIER,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),
                        (8,1,maxx*12*FLOW_DOMAIN_MULTIPLIER,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),(8,2,maxx*12*FLOW_DOMAIN_MULTIPLIER,0,0),(9,1,maxx*12*FLOW_DOMAIN_MULTIPLIER,0,0),(9,2,maxx+xstep,0,0)]
    else:
        file.write("1"+"    "+str(pointer)+"   "+str(maxx+xstep)+"    "+str(maxy)+"   "+"0\n") #Terminates the loop
        domainMatrix = [(2,1,maxx+xstep,maxy,0),(2,2,maxx+xstep,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),(3,1,maxx+xstep,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),(3,2,maxx*12*FLOW_DOMAIN_MULTIPLIER,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),(4,1,maxx*12*FLOW_DOMAIN_MULTIPLIER,maxy*5*FLOW_DOMAIN_MULTIPLIER,0),(4,2,maxx*12*FLOW_DOMAIN_MULTIPLIER,0,0),
                        (5,1,maxx*12*FLOW_DOMAIN_MULTIPLIER,0,0),(5,2,chamberx,0,0),(6,1,chamberx,0,0),(6,2,chamberx,chambery,0)]
        
    for domain_point in range (0, len(domainMatrix)):
        file.write("\n"+str(domainMatrix[domain_point][0])+"    "+str(domainMatrix[domain_point][1])+"   "+str(domainMatrix[domain_point][2])+"    "+str(domainMatrix[domain_point][3])+"   "+str(domainMatrix[domain_point][4])+"\n")

        

    file.close()