import sys

import math

import os

import json

import random
import matplotlib.pylab as plt
import seaborn as sns
import matplotlib
from itertools import permutations

import gurobipy as gp
import numpy as np
import pandas as pd
from mpmath.libmp.backend import xrange
from numpy import sin, cos
from scipy.spatial import distance_matrix
from numpy.random import multinomial

from gurobipy import GRB

def create_data_model():

    """Stores the data for the Pollution Routing Problem."""

    data = {}
    # Initialize a random points matrix with values between 0, 10 (all points in the upper right 0,10 quadrant)
    rnd = np.random
    np.random.seed(7)

    rows = 13 # number of points
    columns = 2  # number of dimensions - 2=2D, 3=3D etc.
    samples = np.empty((rows, columns))
    for i in xrange(0, rows):
        for j in xrange(0, columns):
            samples[i][j] = rnd.randint(0, 300)

    df3 = pd.DataFrame(samples, columns=['x', 'y'])
    distances = pd.DataFrame(distance_matrix(df3[['x', 'y']].values, df3[['x', 'y']].values), index=df3.index,columns = df3.index).values

    data['distance_matrix'] = distances
    data['num_vehicles'] = 5
    data['depot'] = 0
    data['vehicle_capacities'] = [12.0, 12.0, 12.0]
    rn = rows
    rcap = data['vehicle_capacities']
    print (rcap)
    rm = rcap[0] * data['num_vehicles']
    b = rnd.multinomial(rm, np.ones(rn) / rn)
    #data['demands'] = [0,5,12,7,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]
    data['demands'] =b
    data['df3'] = df3

    return data


data = create_data_model()
newdf3 =  data['df3']
distList = data['distance_matrix']

fig, ax = plt.subplots(figsize=(8, 7))
sns.heatmap(distList, ax=ax, cmap='Blues', annot=True, fmt='.0f', cbar=True, cbar_kws={"shrink": .3}, linewidths=.1)
plt.title('distance matrix')
plt.show()


Q = data['demands']

capacities = data['vehicle_capacities']

capacity = capacities[0]

print(Q)

n = len(distList)
dist = {(i, j): distList[i][j] for i in range(n) for j in range(n)}
d = {i: Q[i] for i in range(1,n)}

K = 6

#if K > n:

    #print("npoints must be at least as large as ntrucks")

    #sys.exit(1)

print(f"Total {n-1} customers, {K} trucks with {capacity} capacity")

locations = np.arange(0, 25, 1).tolist()
print(locations[0])


# draw problem state


for i, row in newdf3.iterrows():
    if i == 0:
        plt.scatter(row['x'], row['y'], c='r')
        plt.text(row['x'] + 1, row['y'] + 1, f'Depot - {locations[0]}')
    else:
        plt.scatter(row['x'], row['y'], c='black')
        plt.text(row['x'] + 1, row['y'] + 1, f'{i}({Q[i]}) - {locations[i]}')

plt.xlim([-10, n *50])
plt.ylim([-10, n *50])
plt.title('Kunde: Nummer(Nachfrage)')
plt.show()

################################################################# Create Gurobi Model
m = gp.Model("MultObj")

################################################################# Create variables

x = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='x')
y = m.addVars(d, vtype=GRB.CONTINUOUS, name='y')

################################################################# add truck number as a variable
p = m.addVar(vtype=GRB.INTEGER, obj=1.0, name="p", lb=1, ub=K)

################################################################# Track cumulative demand at each node; cannot exceed capacity
u = m.addVars(n, ub=capacity, name='u')
pairs = [(i, j) for i in range(n) for j in range(1, n) if i != j]
m.addConstrs((u[j] >= u[i] + Q[j] - (1 - x[i, j]) * capacity for (i, j) in pairs), 'demand')

################################################################# Depot cumulative demand is always 0
u[0].LB = 0
u[0].UB = 0
m.update()



# Add the constraints
m.addConstr(gp.quicksum(x[0, j] for j in range(1, n)) == p)  # (9)

inflow = {j: gp.quicksum(x[i, j] for i in range(n) if i != j) for j in range(1, n)}  # (10)
m.addConstrs(inflow[j] == 1 for j in range(1, n))

outflow = {i: gp.quicksum(x[i, j] for j in range(1, n) if j != i) for i in range(1, n)}  # (11)
m.addConstrs(outflow[i] <= 1 for i in range(1, n))

# Inbound and outbound flow is always 1, except for depot
m.addConstrs(x.sum(i, '*') == 1 for i in range(1, n))
m.addConstrs(x.sum('*', i) == 1 for i in range(1, n))



for i in range(1, n):
    for j in range(1, n):
        if i != j:
            m.addConstr(y[i] + d[j] * x[i, j] - capacity * (1 - x[i, j]) <= y[j])

for i in range(1,n):
    m.addConstr(y[i] >= d[i])
    m.addConstr(y[i] <= capacity)
m.update()

input_diesel_leer = 21.3
input_diesel_beladen = 31.7
max_möglich = 12000

#speed = gp.quicksum(((((input_diesel_leer + (input_diesel_beladen - input_diesel_leer) * ((u[i] + Q[j] - capacity * (1 - x[i, j]))/max_möglich)) * 2.629)/100) * dist[i,j]) for i in range(n) for j in range(n))
#speed = gp.quicksum(((((input_diesel_leer + (input_diesel_beladen - input_diesel_leer) * (((u[i] + Q[j] - capacity) * (1 - x[i, j]))/max_möglich)) * 2.629)/100) * dist[i,j]) for i in range(n) for j in range(n))
#speed =  gp.quicksum(dist[i,j] * x[i, j] * Q[j] * 0.078 for i in range(n) for j in range(n))
#totaldist = gp.quicksum(dist[i,j] * x[i,j] for i in range(n) for j in range(n))

#totaldistpart2 = gp.quicksum(0.0981 * Q[i] * dist[i,j] for i in range(n) for j in range(n))

#totaldistpart2 = gp.quicksum(0.0981 * Q[i] * dist[i,j] for i in range(n) for j in range(n))

cd = 0.7
#Air density (kilogram/meter3)
p = 1.2041
#Frontal surface area (meter2)
A = 3.912

#prp_totaldist = gp.quicksum(1.4 *((9.81 * sin(0) + 9.81 * 0.01 * cos(0))* dist[i,j] * 1000 * 26000 * x[i,j]) for i in range(n) for j in range(n))

#prp_totaldist = gp.quicksum(1.4 *(0.0981* dist[i,j] * 1000 * 26000 * x[i,j]) for i in range(n) for j in range(n))
#prp_totaldist = gp.quicksum(1.4 *(9.81 * sin(0) + 9.81 * 0.008 * cos(0))* dist[i,j] * 1000 * 26000 * x[i,j] for i in range(n) for j in range(n))
#prp_loadcost = gp.quicksum((Q[j]*1000) * dist[i,j] * 1000 * x[i,j] for i in range(n) for j in range(n))

#prp_load = gp.quicksum(1.4 * ((9.81 * sin(0) + 9.81 * 0.01 * cos(0)) * dist[i,j] * 1000 * (Q[j])*1000) for i in range(n) for j in range(n))
#prp_totaldist = gp.quicksum(1.4 *(9.81 * sin(0) + 9.81 * 0.008 * cos(0))* dist[i,j] * 1000 * 13000 * x[i,j] for i in range(n) for j in range(n))

#prp_loadcost = gp.quicksum((Q[j]*1000) * dist[i,j] * 1000 * x[i,j] for i in range(n) for j in range(n))

#prp_load = gp.quicksum(1.4 * (9.81 * sin(0) + 9.81 * 0.01 * cos(0)) * dist[i,j] *1000 *  ((u[j]))*1000 for (i, j) in pairs)

#prp_totaldist = gp.quicksum(1.4 *(9.81 * sin(0) + 9.81 * 0.008 * cos(0))* dist[i,j] * 1000 * 26000 * x[i,j] for i in range(n) for j in range(n))
prp_totaldist = gp.quicksum(1.4 *(9.81 * sin(0) + 9.81 * 0.008 * cos(0))* dist[i,j] * 1000 *  26000*x[i, j]  for i in range(n) for j in range(n))

#prp_loadcost = gp.quicksum((Q[j]*1000) * dist[i,j] * 1000 * x[i,j] for i in range(n) for j in range(n))

prp_load = gp.quicksum(1.4 * (9.81 * sin(0) + 9.81 * 0.008 * cos(0)) * dist[i,j] *1000* (38-(u[i] + Q[j]))* 1000*x[i, j]  for i in range(n) for j in range(n))

#(u[i] + Q[j])
#prp_load = gp.quicksum(1.4 * (9.81 * sin(0) + 9.81 * 0.01 * cos(0)) * dist[i,j] * 1000 * x[i, j]* (Q[j])*1000  for (i, j) in pairs)

#prp_load = gp.quicksum(1.4 * (9.81 * sin(0) + 9.81 * 0.008 * cos(0)) * dist[i,j] * 1000  * x[i, j]* (26 - u[j])*1000 for (i, j) in pairs)
#macf = gp.quicksum((Q[j])  * x[i,j] for (i, j) in pairs)

#speed =  gp.quicksum(dist[i,j] * x[i, j] * Q[j] * 0.078 for i in range(n) for j in range(n))


prp_speed = gp.quicksum((1.4) * 0.5 * 1.2041 * 0.7 * 8.2 * dist[i,j] *1000 * (40 * 0.2777778)**2 for i in range(n) for j in range(n))
prp_driver = gp.quicksum(0.0033 * ((dist[i,j]*1000)/(40 * 0.2777778)) for i in range(n) for j in range(n))





td = gp.quicksum(x[i,j] * dist[i,j]  for i in range(n) for j in range(n))

#co2 = gp.quicksum(x[i,j] * ((capacity) - (Q[i])) for i in Q for j in Q)

#kd = gp.quicksum(x[i,j] * (capacity - (Q[i])) for i in Q for j in Q)


#tt = gp.quicksum(0.078 * (capacity - (Q[i]))  * dist[i,j] for i in range(n) for j in range(n) for i in Q)

#prp_loadcost = gp.quicksum(x[i,j]* (7 - Q[i]) * dist[i,j] for (i, j) in pairs)

#m.setObjectiveN(prp_totaldist, 0)
#m.setObjectiveN(-prp_load,1)
#m.setObjectiveN(prp_speed, 2)
#m.setObjectiveN(prp_driver, 2)
#m.setObjectiveN(speed,0)


#co2 = gp.quicksum([x[i,j] * (capacity - (Q[i])) for i in Q for j in Q])

m.setObjective(prp_totaldist + prp_load + prp_speed + prp_driver, GRB.MINIMIZE)
#m.setObjective(td, GRB.MINIMIZE)
#m.setObjective(prp_totaldist + prp_load + prp_speed + prp_driver, GRB.MINIMIZE)
#m.setObjective(speed, GRB.MINIMIZE)
# Optimize multi object model
#m.setObjectiveN(prp_totaldist, 0)
#m.setOjbectiveN(prp_load, 1)

m.update()

m._x = x
m.update()
m.optimize()
# Print optimal routes
gesamt = 0
vals = m.getAttr('X', x)
selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.99)
for i, tup in enumerate(selected.select(0, '*')):
    print("\nRoute for truck {}:\n 0 Load(0)".format(i+1), end='')
    neighbor = tup[1]
    truck_dist = distList[0][neighbor]
    truck_load = Q[neighbor]
    while neighbor:
        print(" -> {} Load({})".format(neighbor, truck_load), end='')
        next_neighbor = selected.select(neighbor, '*')[0][1]
        truck_dist += distList[neighbor][next_neighbor]
        truck_load += Q[next_neighbor]
        neighbor = next_neighbor
    gesamt += truck_dist
    print("gesamt", gesamt)
    print(" -> 0 Load({})".format(truck_load))
    print("Route distance: {}".format(truck_dist))
    print("Route load: {}".format(truck_load))
print("gesamt", gesamt)
print("Route distance: {}".format(truck_dist))
print("\nTotal distance for all routes: {}".format(m.ObjVal))



plt.figure(figsize=(5, 5))

# draw problem state
for i, row in newdf3.iterrows():
    if i == 0:
        plt.scatter(row['x'], row['y'], c='r')
        plt.text(row['x'] + 1, row['y'] + 1, 'depot')
    else:
        plt.scatter(row['x'], row['y'], c='black')
        plt.text(row['x'] + 1, row['y'] + 1,f'{i}({Q[i]}) - {locations[i]}')

plt.xlim([-10, n *50])
plt.ylim([-10, n *50])
plt.title('Kunde: Nummer(Nachfrage)')


# draw optimal route
vals = m.getAttr('X', x)
cmap = matplotlib.cm.get_cmap('Dark2')
selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.99)
for i, tup in enumerate(selected.select(0, '*')):
     # identify the route of each vehicle
    vehicle_route = [selected[i]]

    while vehicle_route[-1][1] != 0:
        for p in selected:
            if p[0] == vehicle_route[-1][1]:
                vehicle_route.append(p)
                # print(vehicle_route)
                break

    # draw for each vehicle
    arrowprops = dict(arrowstyle='->', connectionstyle='arc3', edgecolor=cmap(i))
    for i, j in vehicle_route:
        plt.annotate('', xy=[newdf3.iloc[j]['x'], newdf3.iloc[j]['y']], xytext=[newdf3.iloc[i]['x'], newdf3.iloc[i]['y']],
                         arrowprops=arrowprops)
plt.show()
# Print optimal routes
# km/h in m/s Umrechnung
pay = 10
fuel_price = 1.4
v = 40 * 0.2777778
vkmh = 40
convfactor = 8.8
empty_vehicle = 7.5
beta = 0.5 * 0.7 * 5 * 1.2041
energy = 0
gesamt = 0
fuel = 0
gfuel = 0
#Fuel-to-air mass ratio
fam = 1
#Heating value of a typical diesel fuel (kilojoule/gram
heat = 44
#Conversion factor (gram/second to liter/second)
conv = 737

# lambda - first part of the fuel consumption function
la = fam/(heat * conv)

#Engine friction factor (kilojoule/rev/liter)
k = 0.2
#Engine speed (rev/second)
N = 33
#Engine displacement (liters)
V = 5

kNV = k*N*V

# w  as empty vehicle weight in kilogram
w = 6350

grav = 9.81

#y1 = la * (kNV * 100000/(x*0.277778))
# speed module
#y2 = la * ((1/(1000*0.4*0.9)) * (0.5*0.7*4*1.2 * 100000 * ((x*0.277778)**2)))
# weight module
#y3 = la * ((1/(1000*0.4*0.9)) * (grav * sin(0) + grav * 0.01 * cos(0)) * (6350 + x)*100000)

#y4 = y1 + y2 + y3

vals = m.getAttr('X', x)
selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.99)
if empty_vehicle <= 10:
    for i, tup in enumerate(selected.select(0, '*')):
        print("\nRoute for truck {}:\n 0 Load(6)".format(i+1), end='')
        neighbor = tup[1]
        truck_dist = distList[0][neighbor]
        #truck_load = Q[neighbor]
        truck_typ_load = capacity
        print(capacity)
        #if truck_typ_load <= 4.0:

        truck_load = capacity + 17
        #truck_load = capacity + 3

        truck_loadminus = Q[neighbor]
        #load-induced energy requirements in Joule // Energie durch 3,611e+6 dividieren, um ein ungefähres Ergebnis zu erhalten
        load_induced = distList[0][neighbor] * 1000 * 0.0981 * truck_load * 1000
        #speed-induced energy requirements in Joule // Energie durch 3,611e+6 dividieren, um ein ungefähres Ergebnis zu erhalten
        speed_induced = distList[0][neighbor] * 1000 * v*v * beta
        weight_fuel_rate = (la * (kNV * 100000/(v)) + la * ((1/(1000*0.4*0.9)) * (0.5*0.7*4*1.2 * 100000 * ((v)**2))) + la * ((1/(1000*0.4*0.9)) * (grav * sin(0) + grav * 0.01 * cos(0)) * (truck_load*1000)*100000))/100 * distList[0][neighbor]
        print("Liter nach Gewicht", weight_fuel_rate)
        route_time = (distList[0][neighbor])/vkmh
        KwRoute = (load_induced + speed_induced) * route_time
        frate = (KwRoute/3600000)/convfactor
        averageco2 = (((input_diesel_leer + (input_diesel_beladen - input_diesel_leer) * (truck_load*1000 / max_möglich)) * 2.629) / 100) * (distList[0][neighbor])
        averagefuel = ((input_diesel_leer + (input_diesel_beladen - input_diesel_leer) * (truck_load*1000 / max_möglich))/ 100) * (distList[0][neighbor])
        print("teeeeeest fuel", KwRoute,route_time, frate)
        print("erste route", load_induced, "aus", distList[0][neighbor],truck_load)
        print("averagefuel", averageco2, averagefuel)
        print("speed_induced route", speed_induced, "aus", distList[0][neighbor], truck_load)
        while neighbor:
            truck_load -= Q[neighbor]
            print(" -> {} Load({})".format(neighbor, truck_load, ), end='')
            next_neighbor = selected.select(neighbor, '*')[0][1]
            tr_dist = distList[neighbor][next_neighbor]
            print("While Schleife Route distance: {} zu multiplizierendes Gewicht {}".format(truck_dist,truck_load))
            truck_dist += distList[neighbor][next_neighbor]
            nextroute = distList[neighbor][next_neighbor] * 1000 * 0.0981 * truck_load * 1000
            nextroutes_speed_induced = distList[neighbor][next_neighbor] * 1000 * v*v * beta
            next_route_time = (distList[neighbor][next_neighbor]) / vkmh
            next_averagefuel = ((input_diesel_leer + (input_diesel_beladen - input_diesel_leer) * (truck_load *1000/ max_möglich)) / 100) * (distList[neighbor][next_neighbor])
            next_weight_fuel_rate = (la * (kNV * 100000 / (v)) + la * (
                        (1 / (1000 * 0.4 * 0.9)) * (0.5 * 0.7 * 4 * 1.2 * 100000 * ((v) ** 2))) + la * (
                                            (1 / (1000 * 0.4 * 0.9)) * (grav * sin(0) + grav * 0.01 * cos(0)) * (
                                                 truck_load * 1000) * 100000)) / 100 * (distList[neighbor][next_neighbor])
            weight_fuel_rate += next_weight_fuel_rate
            print("Sprit nach Formel Bectas",  weight_fuel_rate)
            print("time",next_route_time )
            next_KwRoute = (nextroute + nextroutes_speed_induced) * next_route_time
            print("next_KwRoute", next_KwRoute)
            next_frate = (next_KwRoute / 3600000) / 8.8
            print("fuel consumption next_frate", next_frate)
            averagefuel += next_averagefuel
            load_induced += nextroute
            speed_induced += nextroutes_speed_induced
            frate += next_frate
            route_time += next_route_time
            print("fuel consumption", frate, route_time)
            print("nextroute", nextroute)
            print("nextroutes_speed_induced",nextroutes_speed_induced)
            print("load_induced", load_induced)
            print("speed_induced", speed_induced)
            neighbor = next_neighbor

        gesamt += truck_dist
        energy = load_induced + speed_induced
        fuel += averagefuel
        gfuel += weight_fuel_rate
        #print("While Schleife Route distance: {} zu multiplizierendes Gewicht {} {}".format(truck_dist, truck_load))
        print(" -> 0 Load({})".format(truck_load))
        print("Route distance: {}".format(truck_dist))
        print("Route load: {}".format(truck_load))
    print("total fuel consumpted in liter",averagefuel)
    print("fuel consumption", frate,"Spritkosten", frate * fuel_price,  "time spent",route_time, "Pay", route_time * pay, "Co2 Kosten",((frate * 2.32)/1000)*27 )
    print("load_induced", load_induced)
    print("speed_induced", speed_induced)
    print("Gesamtkosten", energy)
    #print("While Schleife Route distance: {} zu multiplizierendes Gewicht {} {}".format(truck_dist, truck_load, tr_dist))
    print("\nobjectvalue: {}".format(m.ObjVal))
    print("gesamt", gesamt,"Zeit", gesamt / 40 , "Fahrerkosten", (gesamt / 40)*12,
          "\nfuel", fuel, "FuelConsumption nach Bectas", gfuel, "Kosten Sprit",gfuel*1.40,
          "\nCO2 Wert", (gfuel * 2.629)/1000, "CO2 Kosten", ((gfuel * 2.629)/1000)*25,
          "\nGesamtkosten", (gesamt / 40)*12 + gfuel*1.40 + ((gfuel * 2.629)/1000)*25)



#m.write("out.mst")
#m.write("out.sol")
#m.write("out.lp")

