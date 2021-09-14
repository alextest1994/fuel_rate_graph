from math import sin, cos

import matplotlib.pyplot as plt
import numpy as np
from numpy import interp
from scipy.optimize import fmin
import math
from scipy.optimize import minimize
# 100 linearly spaced numbers
x = np.linspace(0,10000,10000)


#Fuel-to-air mass ratio
fam = 1
#Heating value of a typical diesel fuel (kilojoule/gram
heat = 44
#Conversion factor (gram/second to liter/second)
conv = 737

# lambda - first part of the fuel consumption function
la = fam/(heat * conv)

print(la)
#kNV - will only be significant for low speed levels
#Engine friction factor (kilojoule/rev/liter)
k = 0.2
#Engine speed (rev/second)
N = 33
#Engine displacement (liters)
V = 5

kNV = k*N*V

# w  as empty vehicle weight in kilogram
w = 6350
# 0 - 3650 Payload
q = 3650
# constant gamma

# xinst
xinst = 0

rogr = 0
#ntf Vehicle drive train efficiency = 0.4
#n Efficiency parameter for diesel engines = 0.9
gamma = 1/(1000 * 0.4 * 0.9)

nef = 0.9
#gamma = 3
#Coefficient of rolling resistance
cr = 0.01
#vehicle-arc specific constant
alpha = 9.81 * cr

grav = 9.81
#vehicle load
f = 1500


#Coefficient of aerodynamic drag
cd = 0.7
#Air density (kilogram/meter3)
p = 1.2041
#Frontal surface area (meter2)
A = 3.912
#vehicle specific constant beta
beta = 0.5 * cd * p * A
ntf = 0.4
pacc = 0
# Kosten pro Liter
cf = 1

termnenner = 1000 * heat * conv * ntf * nef
# the function
input_diesel_leer = 21.3
input_diesel_beladen = 31.7
beladung = 4000
max_möglich = 4000
#y = la*(0.2 * 33 * 5 + w * gamma * alpha *x + gamma * alpha * x * f + beta * gamma *x*x*x )*(100/(x))
#y = (la* (33 + w * gamma * 0.0981 * x + gamma * 0.0981 * x * f + 2.1 * gamma * (x)**3)) *10000/(x)
#y = 0.5*0.7*5*1.2041*x**3 + 7000 * 9.81 * 0.01 * x


#obere Zeile ohne v multiplikation
term1 = (((w+ q)*(xinst + grav * 0* rogr + grav * cr * 1) + 0.5 * cd * p * A * ((x*0.277778)**2))*(x*0.277778))

#y = fam/(heat * conv) * (kNV + (1/nef) * ((((w+ q)*(xinst + grav * 0* rogr + grav * cr * 1) + 0.5 * cd * p * A * x*x)*x) / (1000 * ntf) + pacc))
c0 = (cf * w * 1 * (xinst + grav * 0* rogr + grav * cr * 1)) / (termnenner)
c1 = c0 / w
c2 = (cf * 1 * kNV) / (heat * conv)
c3 = (1 * beta)/(termnenner)
# v = optimale
vopt = (c2 / (2 * c3)) ** (1/3)

print(vopt)
# setting the axes at the centre
# fuel consumption in euro

d = 1000

#y = (c0 + c2/x + c3 * ((x*0.277778)**2)) * 100000 + c1 * 100000

x = np.linspace(0,10000,10000)

#y1 = la*(kNV * 100000/(x))
# engine modul
y1 = la * (kNV * 100000/(x*0.277778))
# speed module
y2 = la * ((1/(1000*0.4*0.9)) * (0.5*0.7*4*1.2 * 100000 * ((x*0.277778)**2)))
# weight module
y3 = la * ((1/(1000*0.4*0.9)) * (grav * sin(0) + grav * 0.01 * cos(0)) * (x)*100000)

y4 = y1 + y2 + y3
y5 = y3
y = la * (kNV * 100000/(x*0.277778)) + la * ((1/(1000*0.4*0.9)) * 0.5*0.7*4*1.2 * 100000 * ((x*0.277778)**2)) + la * ((1/(1000*0.4*0.9)) * (grav * sin(0) + grav * 0.01 * cos(0)) * (6350 + x)*100000)
np.all(np.diff(x) > 0)
np.all(np.diff(y) > 0)

x_value = interp(15, y, x)

print("test", x_value)
#y = la*(gamma * beta * 100000 * x**2)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(5, 10000)
ax.set_ylim([0, 30])
# plot the function
plt.plot(x, y5, 'b', label="Gewicht")
#plt.plot(x, y2, 'r', label="Geschwindigkeitskomponente")
#plt.plot(x, y3, 'g', label="Gewichtskomponente")
#plt.plot(x, y4, 'purple', label="Kraftstoffverbrauch")
#plt.plot(x, y, 'p', label="test")

print(np.interp(15, y4, x))

plt.title("Kraftstoffverbrauch",fontsize=15)
plt.xlabel("Gewicht in Kilogramm",fontsize=13)
plt.ylabel("F(Liter/100km)",fontsize=13)
plt.legend()
plt.show()
# show the plot
plt.show()


print(y)