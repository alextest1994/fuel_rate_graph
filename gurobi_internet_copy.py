import sys
#test
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
from scipy.spatial import distance_matrix


from gurobipy import GRB

def make_edges(ll):

    tuples = []

    for ii in range(1, len(ll)):

        tuples.append((ll[ii-1], ll[ii]))

    return tuples

# Callback - use lazy constraints to eliminate sub-tours

def subtourelim(model, where):

    if where == GRB.Callback.MIPSOL:

        # make a list of edges selected in the solution

        vals = model.cbGetSolution(model._x)

        selected = gp.tuplelist((i, j) for i, j in model._x.keys()

                                if vals[i, j] > 0.5)

        # add the capacity constraints

        for i, tup in enumerate(selected.select(0, '*')):

            capacity_k = 0

            nodes_k = [0]

            neighbor = tup[1]

            while neighbor:

                capacity_k += Q[neighbor]

                nodes_k.append(neighbor)

                neighbor = selected.select(neighbor, '*')[0][1]

            if capacity_k > capacity:

                model.cbLazy(gp.quicksum(Q[j]*model._x[i,j] for i, j in make_edges(nodes_k)) <= capacity)

        # find the shortest cycle in the selected edge list

        tour = subtour(selected)

        if len(tour) < n:

            # add subtour elimination constr. for every pair of cities in tour

            model.cbLazy(gp.quicksum(model._x[i, j]

                                     for i, j in permutations(tour, 2))

                         <= len(tour)-1)

# Given a tuplelist of edges, find the shortest subtour not containing depot

def subtour(edges):

    unvisited = list(range(1, n))

    cycle = range(n+1)  # initial length has 1 more city

    # First, remove all nodes connected to depot

    depot_connected = [j for i, j in edges.select(0, '*') if j != 0]

    # print('Depot_connected:', depot_connected)

    while depot_connected:

        current = depot_connected.pop()

        # print('Current:', current)

        # print('Unvisited:', unvisited)

        unvisited.remove(current)

        neighbors = [j for i, j in edges.select(current, '*')

                     if j in unvisited and j != 0]

        depot_connected += neighbors

    # Now, find subtour using tsp.py code

    while unvisited:

        thiscycle = []

        neighbors = unvisited

        while neighbors:

            current = neighbors[0]

            thiscycle.append(current)

            unvisited.remove(current)

            neighbors = [j for i, j in edges.select(current, '*')

                         if j in unvisited]

        if len(cycle) > len(thiscycle):

            cycle = thiscycle

    return cycle

def create_data_model():

    """Stores the data for the Pollution Routing Problem."""

    data = {}
    df = pd.DataFrame(np.array([[55, 55], [35, 60], [10, 60], [50, 70], [75, 90], [100, 40], [70, 10], [40, 20], [40, 40], [20, 10], [5, 15]]), columns=['x', 'y'])
    distances = pd.DataFrame(distance_matrix(df[['x', 'y']].values, df[['x', 'y']].values), index=df.index,
                             columns=df.index).values

    data['distance_matrix'] = distances
    data['num_vehicles'] = 3
    data['depot'] = 0
    data['demands'] = [0.0, 3.0, 3.0, 6.0, 3.0, 3.0, 6.0, 4.0, 2.0, 2.5, 3.5]
    data['vehicle_capacities'] = [12.0, 12.0, 12.0, 12.0]

    return data


data = create_data_model()

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

K = data['num_vehicles']

if K > n:

    print("npoints must be at least as large as ntrucks")

    sys.exit(1)

print(f"Total {n-1} customers, {K} trucks with {capacity} capacity")

locations = ['Jena','Weimar', 'Erfurt','Apolda','Naumburg','Gera','Neustadt','test1','test2','test3','test4']
print(locations[0])


# draw problem state
df = pd.DataFrame(np.array([[55, 55, 0], [35, 60, 3],
                            [10, 60, 3], [50, 70, 6],
                            [75, 90, 3], [100, 40, 3],
                            [70, 10, 6], [40, 20, 4],[40, 40, 2],[20, 10, 2.5],[5, 15, 3.5]]),
                            columns=['x', 'y', 'demand'])
for i, row in df.iterrows():
    if i == 0:
        plt.scatter(row['x'], row['y'], c='r')
        plt.text(row['x'] + 1, row['y'] + 1, f'Depot - {locations[0]}')
    else:
        plt.scatter(row['x'], row['y'], c='black')
        plt.text(row['x'] + 1, row['y'] + 1, f'{i}({Q[i]}) - {locations[i]}')

plt.xlim([-10, 110])
plt.ylim([-10, 110])
plt.title('points: id(demand)')
plt.show()

m = gp.Model()

# Create variables

x = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='e')

# Inbound and outbound flow is always 1, except for depot

m.addConstrs(x.sum(i, '*') == 1 for i in range(1, n))

m.addConstrs(x.sum('*', i) == 1 for i in range(1, n))

# all customers are served

AllNodesVisited={i: gp.quicksum(x[j,i] for j in range(n) if i!=j) for i in range(1,n)} # (2.2)

m.addConstrs(AllNodesVisited[i]==1 for i in range(1,n))

# Depot has inbound and outbound flow equal to number of trucks

m.addConstr(x.sum(0, '*') <= K)

m.addConstr(x.sum('*', 0) <= K)

# Track cumulative demand at each node; cannot exceed capacity
u = m.addVars(n, ub=capacity, name='u')
pairs = [(i, j) for i in range(n) for j in range(1, n) if i != j]
m.addConstrs((u[j] >= u[i] + Q[j] - capacity * (1 - x[i, j])
             for (i, j) in pairs), 'demand')

# Depot cumulative demand is always 0
u[0].LB = 0
u[0].UB = 0

input_diesel_leer = 21.3
input_diesel_beladen = 31.7
max_möglich = 12000

#speed = gp.quicksum(((((input_diesel_leer + (input_diesel_beladen - input_diesel_leer) * (((u[i] + Q[j]))/max_möglich)) * 2.629)/100) * dist[i,j]) * (x[i, j])for i in range(n) for j in range(n))
#speed = gp.quicksum(((((input_diesel_leer + (input_diesel_beladen - input_diesel_leer) * ((u[i] + Q[j] - capacity * (1 - x[i, j]))/max_möglich)) * 2.629)/100) * dist[i,j]) for i in range(n) for j in range(n))
speed =  gp.quicksum(dist[i,j] * x[i, j] * Q[j] * 0.078 for i in range(n) for j in range(n))
#totaldist = gp.quicksum(dist[i,j] * x[i,j] for i in range(n) for j in range(n))

#totaldistpart2 = gp.quicksum(0.0981 * Q[i] * dist[i,j] for i in range(n) for j in range(n))

#totaldistpart2 = gp.quicksum(0.0981 * Q[i] * dist[i,j] for i in range(n) for j in range(n))

#prp_totaldist = gp.quicksum(dist[i,j] * x[i,j]   for i in range(n) for j in range(n))

#prp_loadcost = gp.quicksum(x[i,j]* (19 - Q[i]) * dist[i,j] for (i, j) in pairs)

#prp_speed = gp.quicksum((1.4 +  5) * dist[i,j] for i in range(n) for j in range(n))

#prp_driver = gp.quicksum(15 * (dist[i,j]/70) for i in range(n) for j in range(n))




#td = gp.quicksum(x[i,j] * dist[i,j] for i in range(n) for j in range(n))

#co2 = gp.quicksum([x[i,j] * ((capacity) - (Q[i])) for i in Q for j in Q])

#kd = gp.quicksum(x[i,j] * (capacity - (Q[i])) for i in Q for j in Q)


#tt = gp.quicksum(0.078 * (capacity - (Q[i]))  * dist[i,j] for i in range(n) for j in range(n) for i in Q)


#m.setObjectiveN(totaldist, 0)
#m.setObjectiveN(prp_loadcost, 1)


#co2 = gp.quicksum([x[i,j] * (capacity - (Q[i])) for i in Q for j in Q])

m.setObjective(speed, GRB.MINIMIZE)



#m.ModelSense = GRB.MINIMIZE
# Optimize model

m._x = x

m.Params.LazyConstraints = 1

m.optimize(subtourelim)

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
for i, row in df.iterrows():
    if i == 0:
        plt.scatter(row['x'], row['y'], c='r')
        plt.text(row['x'] + 1, row['y'] + 1, 'depot')
    else:
        plt.scatter(row['x'], row['y'], c='black')
        demand = row['demand']
        plt.text(row['x'] + 1, row['y'] + 1, f'{i}({demand})')

plt.xlim([-10, 110])
plt.ylim([-10, 110])
plt.title('points: id(demand)')


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
        plt.annotate('', xy=[df.iloc[j]['x'], df.iloc[j]['y']], xytext=[df.iloc[i]['x'], df.iloc[i]['y']],
                         arrowprops=arrowprops)
plt.show()
# Print optimal routes

total_route_length = 0
total_route_costs = 0
total_cost = 0

vals = m.getAttr('X', x)
selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.99)
truck_load_route = capacity
for i, tup in enumerate(selected.select(0, '*')):
    truck_full_loaded = (capacity + 6)
    print("\nRoute for truck {}:\n 0 Load(0) Route_Load {}".format(i+1, truck_full_loaded), end='')
    neighbor = tup[1]
    truck_dist = distList[0][neighbor]
    truck_load = Q[neighbor]
    c = 0
    route_costs = distList[0][neighbor] * 18 * 0.078
    print("RTC", route_costs)
    tc = 0
    while neighbor:
        truck_full_loaded -= Q[neighbor]
        route_costs += truck_full_loaded * truck_dist * 0.078
        print("\nKosten mit 0.078: {}".format(total_cost))
        print(" -> {} Load({}) Route_Load {}".format(neighbor, truck_load, truck_full_loaded), end='')
        next_neighbor = selected.select(neighbor, '*')[0][1]
        #route_costs += distList[neighbor][next_neighbor]

        truck_dist += distList[neighbor][next_neighbor]
        c += distList[neighbor][next_neighbor] * truck_full_loaded
        tc += c * 0.078
        print("\nccc", c, tc)
        truck_load += Q[next_neighbor]
        neighbor = next_neighbor

        #total_cost = total_cost + route_costs
    total_route_length += truck_dist
    #total_cost = distList[neighbor][0]
    total_cost += route_costs
    print("Kosten: {}".format(truck_dist))
    print(" -> 0 Load({})".format(truck_load))
    print("Route distance: {}".format(truck_dist))
    print("Route load: {}".format(truck_load))
    print("Route load current location: {}".format(truck_load_route))
    print("Route CCCCCCCCCCCCCCCCC: {}".format(c))


print("Total costsdist different calculation: {}".format(total_cost))
print("Total Distance different calculation: {}".format(total_route_length))
print("\nTotal distance for all routes: {}".format(m.ObjVal))
print("Route TTTTTCCCCCCCCCCCCCCCCC: {}".format(c))

m.write("out.mst")
m.write("out.sol")
m.write("out.lp")