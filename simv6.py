"""
    V6
    -Runge-Kutta for heat transfer
    First update: 03/08/2024
    Last update: 03/08/2024
"""

import math
from math import pi
import matplotlib.pyplot as plt   #used for plotting different things 
from os import remove
from os.path import exists

fnames = ["Material Variation Table.txt", "CombDyn Table.txt", "performance.txt"]

for f in range(0, len(fnames)):
    if (exists(fnames[f])):  #Removes old version of this table
        remove(fnames[f])

inch = 0.0254

markers = ['o', 'v', '^', 's', '8']
linestyles = ['-',  ':']

def linstl(i):
    j=0
    while i> len(markers)-1:
        i -= len(markers)
        j+=1
    return linestyles[j]
def mrkr(i):
    while i>len(markers)-1:
        i -= len(markers)
    return markers[i]



"""*****************************************************"""
"""SETTINGS"""
COMBSim = True
PressAtm = True             #Exit pressure is atmospheric (101325 Pa) 
VARPressure = True          #Requires COMBSim, pressure step size and amount
VARMaterials = True
TWRSim = False              #Requires rocket mass
CALCMass = True             #Requires VARMaterials
ROCKETSim = True
ROCKETISPVar = True         #Requires VARPressure

TEMPSim = True             #Requires VARMaterials

graphs = [True , False, False, False, False, False, False, False, False, True]
        # x,     v,     a,     m,     gl,    dl,    tl,    re,    ma,    T_s
"""*****************************************************"""

chamberpressure = 5000000   #test pressure to be used if it is not to vary.
pressincrease = 1000000     #Pressure step (i.e. 50k 100k 150k 200k -> step is 50k)
presssteps = 6              #Number of pressure steps (above example: 4)
desiredThrust = 600         #Thrust to be used if it is not to vary.
desiredFuelRatio = 0.6      #The ratio between fuel mass and total mass.
maxHeight = 0.7               #Maximum height if mass ratio is not met.

desiredTWR = 3.0            #This is used in the ROCKETSim module to calculate initial thrust if all other calculations are not realized.

g0 = 9.81           #constant used for calculating specific impulse
Cd = 0.5            #drag coefficient - not yet varying

m_payload = 1.0
m_nozzle = 2.0

FS = 2              #safety factor for casing thickness
bltz = 5.67*(10**-8)    #W/m^2K^4

"""*****************************************************"""
"""DIMENSIONS"""
d_casing = 4*inch       #outer casing diameter is 4 in
d_fuel = .09             #fuel cylinder diameter is 87.73 mm (autoCAD)
l_payload = [1*d_casing, 2*d_casing, 3*d_casing]    #possible values
l_end = 0.082                                       #from autoCAD
E_s = 0.7                   #weld efficiency
t_insulation = 0.002                #2.0 mm insulation
chosenMaterial = 1      #For casing calculation. See materials[0].
L0 = 0.1                            #Initial length between bottom of fuel cylinder and nozzle [m] (guess)
"""*****************************************************"""

"""*****************************************************"""
"""FUNCTIONS"""

def Re(rho, v, L, mu):      #Calculates Reynolds' number
    return abs(rho*v*L/mu)
def Ma(v, c):               #Calculates Mach number for any given speed
    return abs(v/c)
def linMass(d_o, d_i, rho):
    return(pi/4*rho*(d_o**2 - d_i**2))
def casingTCalc(P, FS, d, sigma):
    r = d/2
    return(FS*P*r/(E_s*sigma*10**6 -0.6*P))
def linintrp (xaxis, yaxis, xinterp):       #Linear interpolation
    if xinterp <= xaxis[0]:                     #if below lowest x, return lowest y
        return yaxis[0]    
    if xinterp > xaxis[-1]:
        return yaxis[-1]
    for x in range(0, len(yaxis)):              #if above highest x, return highest y
        if (xinterp >= xaxis[x] and xinterp <xaxis[x+1]):   #else interpolate
            x1 = xaxis[x]       #x1 is element x[i]
            y1 = yaxis[x]       #y1 is element y[i]
            x2 = xaxis[x+1]     #x2 is element x[i+1]
            y2 = yaxis[x+1]     #y2 is element y[i+1]
            return (y1 + (y2-y1)/(x2-x1)*(xinterp-x1))
def writetable(filename, colHeading, rowHeading, rowUnit, rowScalef, data, dataScalef, spacing, colUnit = None, title = None):         
    f = open(filename, "a")
    if title is not None:
        f.write(" ".ljust(spacing))
        f.write(title + "\n")
    f.write("".rjust(spacing))
    for i in range(0, len(colHeading)):             #Writes all column titles (e.g. Pressure, Area, Diameter...)
        f.write(colHeading[i].ljust(spacing))
    f.write("\n" + rowUnit.ljust(spacing))                         #Writes the unit of the rows
    if colUnit is not None:                         #Checks if columns have units or not
        for i in range(0, len(colUnit)):            #Writes all column units (if there are any)
            f.write(colUnit[i].ljust(spacing))
    for i in range(0, len(rowHeading)):                      #Writes pressure in the first column
        f.write("\n" + str(rowScalef*rowHeading[i]).ljust(spacing))
        for j in range(0, len(colHeading)):
            f.write("{:.3f}".format(dataScalef*data[i*len(colHeading) + j]).ljust(spacing))   #Writes casing thickness for each material and pressure
    f.write("\n\n")                                   #In order to separate from the next table
    f.close()
def combcalc(T_1, k, R, P_1, P_2, F, a, n, rho_fuel):
    T_t = 2*T_1/(k+1)
    V_t = (2*k/(k+1)*R*T_t)**0.5 #Throat velocity. Mach 1. Uses throat temperature.
    
    P_t = P_1*(2/(k+1))**(k/(k-1))
    vol_1 = R*T_t/P_t
    vol_t = (vol_1*(((k+1)/2)**(1/(k-1))))
    V_2 = (((2*k/(k-1)*R*T_1*(1-(P_2/P_1)**((k-1)/k))))**0.5)
    m_dot = (F/V_2)
    A_t = (m_dot*vol_t/V_t)
    AR = (((k+1)/2)**(1/(k-1))*(P_2/P_1)**(1/k)*((k+1)/(k-1)*(1-(P_2/P_1)**((k-1)/k)))**0.5)
    A_2 = (A_t/AR)    
    d_t = ((4/pi*A_t)**0.5)    
    d_2 = ((4/pi*A_2)**0.5)    
    ISP = (V_2/g0)    
    Thr = (V_t*A_t*V_2/vol_t)
    burn_rate = (a*(P_1/1000000)**n)/1000
    A_comb = (m_dot/(rho_fuel*burn_rate))
    d_comb = ((4/pi*A_comb)**0.5)
    
    return [T_t, V_t, P_t, vol_1, vol_t, V_2, m_dot, A_t, AR, A_2, d_t, d_2, ISP, Thr, burn_rate, A_comb, d_comb]
def RKT(h_int, h_ext, T_int, T_w, T_ext, mline, cpmat, t_step):
    def qconv(h, T1, T2):
        return(h*(T1-T2))
    def qrad(T1, T2):
        return(bltz*((T1)**4-T2**4))
    def qtot(T_input):
        q1 = qconv(h_int, T_int, T_input)
        q2 = qrad(T_int, T_input)
        q3 = qconv(h_ext, T_input, T_ext)
        q4 = qrad(T_input, T_ext)
        return(q1 + q2 - q3 - q4)
    
    #/(mline*cpmat)*t_step
    
    k1 = qtot(T_w)/(mline*cpmat)*t_step
    k2 = qtot(T_w + k1*t_step/2)/(mline*cpmat)*t_step/2
    k3 = qtot(T_w + k2*t_step/2)/(mline*cpmat)*t_step/2
    k4 = qtot(T_w + k3*t_step)/(mline*cpmat)*t_step
    
    return (T_w + t_step/6*(k1 + 2*k2 + 2*k3 + k4)/(mline*cpmat))
    
    
    
def trajSim(mf, me, ISP, m_dot, t_step, d_e, d_i, d_insulation, rho_mat, Cp_mat, k_insulation, k_mat, r):
    m = [mf]                    #initial mass
    
    W = [g0*m[0]]               #initial weight
    T = [ISP*m_dot*g0]   #initial thrust
    D = [0.0]                   #initial drag
    FR = [T[0]-W[0]-D[0]]
    Re_t = [0.0]
    Ma_t = [0.0]
    deltaV = [ISP*g0*math.log(m[0])/me]
    
    time = [0.0]        #initial time
    
    a = [FR[0]/m[0]]       #initial acceleration
    v = [0.0]           #initial speed
    x = [0.0]           #initial position
    g_loss = [0.0]
    D_loss = [0.0]
    tot_loss = [0.0]
    
    if TEMPSim:
        T_int = T_1
        T_ext = [linintrp(h, T_atm, x[0])]
        T_s = [T_ext[0]]       #Initial temperature of steel, K
        mline = pi/4*(d_e**2 - d_i**2)*rho_mat
        A_1 = pi/4*((d_i-2*t_insulation)**2)
        L = [L0]
        
        hconv = [3.075*Cp_f*((m_dot/A_1)**0.8)/(d_insulation)*0.2*(1+(A_1/L[0])**0.7)]
        Rh = [1/(pi/4*hconv[0]*d_insulation)]
        Rk = math.log(d_i/d_insulation)/(2*pi*k_insulation)
        
        Nu = 0.140*(Re_t[0])**0.686     #source: Wiberg, Lior (2004)
        h_ext = Nu*k_mat/d_e
        
        Req = [Rh[0] + Rk]
        
    t=0
    while x[t]>=0:
        time.append(time[t]+tstep)
        if (m[t]<=me):  #i.e. if the rocket ran out of fuel
            T.append(0)
            m.append(m[t])
            g_loss.append(g_loss[t])
            if TEMPSim:
                T_int = 298         #when fuel burns out, temperature instantly goes to 25 C
        else:
            T.append(T[0])           #if i ever do change the thrust calculation, this may be handy
            m.append(m[t]- (m_dot*tstep))
            g_loss.append(g_loss[t]+W[t]/m[t]*tstep)

        W.append(m[t]*linintrp(h, g, x[t]))
        
        if v[t]>=0:     #if going up
            D.append(linintrp(h, rho_atm, x[t])*A_transv*Cd*v[t]**2/2)     #drag is in -x
            D_loss.append(D_loss[t]+D[t]*tstep/m[t])
        else:
            D.append(-linintrp(h, rho_atm, x[t])*A_transv*Cd*v[t]**2/2)    #drag is in +x
            D_loss.append(D_loss[t])
        if TEMPSim:
            T_ext.append(linintrp(h, T_atm, x[t]))
            L.append(L[t] + r*t_step)
            hconv.append(3.075*Cp_f*((m_dot/A_1)**0.8)/(d_insulation)*0.2*(1+(A_1/L[t])**0.7))
            Rh.append(1/(pi/4*hconv[t]*d_insulation))
            
            Req.append(Rh[t] + Rk)
            
            Nu = 0.140*(Re_t[t])**0.686     #source: Wiberg, Lior (2004)
            h_ext = Nu*k_mat/d_e
            if h_ext > 50000:
                h_ext = 50000
                        
            T_s.append(RKT(1/Req[t], h_ext, T_int, T_s[t], T_ext[t], mline, Cp_mat, t_step))
            
        FR.append(T[t]-W[t]-D[t])
        a.append(FR[t]/m[t])
        v.append(v[t]+a[t]*tstep)
        x.append(x[t]+v[t]*tstep)
        Re_t.append(Re(linintrp(h, rho_atm, x[t]), v[t], d_casing, linintrp(h, mu, x[t])))
        Ma_t.append(Ma(v[t], linintrp(h, SS_h, x[t])))
        
        tot_loss.append(g_loss[t] + D_loss[t])
        deltaV.append(ISP*g0*math.log(m[t])/me)
        t +=1
    if TEMPSim:
        return [time, x, v, a, m, g_loss, D_loss, tot_loss, Re_t, Ma_t, T_s]
    else:
        return [time, x, v, a, m, g_loss, D_loss, tot_loss, Re_t, Ma_t]

"""*****************************************************"""
"""*****************************************************"""
"""PROPERTIES"""

rho_fuel = 1800 #Fuel density (Nakka)
Cp_f = 1.8134   #Specific heat of sugar candy [kJ/kgK] (Nakka)
k_paper = 0.18  #Thermal conductivity (W/mK)
P_atm = 101325
h = []          #Height (from sea level)
T_atm = []      #Atmospheric temperature (still unused)
g = []          #Gravitational acceleration (could be calculated instead)
P = []          #Atmospheric pressure
rho_atm = []    #Atmospheric density
mu = []         #Dynamic viscosity
SS_h = []       #Speed of sound

f = open("properties.txt", "r")         #Reads the properties file
f.readline()
f.readline()                            #skips the 2 lines of heading
for line in f:                      
    spl = line.split("\t")              #Gets an array of all properties in that line
    h.append(float(spl[0]))             #Gets height
    T_atm.append(float(spl[1]) + 273)   #Gets atmospheric temperature (K)
    g.append(float(spl[2]))             #Gets local grav. acceleration
    P.append(float(spl[3])*10**4)       #Gets local pressure
    rho_atm.append(float(spl[4]))       #Gets local atm. density
    mu.append(float(spl[5])*10**-5)     #Gets local dyn. viscosity
    SS_h.append(float(spl[6]))          #Gets local speed of sound
f.close()

"""*****************************************************"""
"""COMBUSTION PRESSURE VARIATION SETUP"""
if VARPressure:
    P_1 = [i*pressincrease for i in range(1,presssteps+1)]    #Combustion chamber pressure (skips 0 Pa pressure)
else:
    P_1 = [chamberpressure]

n_P = len(P_1)
"""*****************************************************"""

"""****************************************************"""
"""NOZZLE CALCULATIONS"""
R_u = 8134          #Universal gas constant     
MM = 41.98          #Fuel molar mass            (Nakka)
R = R_u/MM          #Specific gas constant
k = 1.1331          #Specific heat ratio        (Nakka)  
T_1 = 1720          #Chamber temperature        (Nakka)
c_star = (R*T_1/(k*(2/(k+1))**((k+1)/(k-1))))**0.5   #characteristic exhaust velocity
a = 8.2635            #Burn rate coefficient (mm/s MPA^n)
n = 0.319           #burn rate coefficient (adimensional) => r = aP^n
A_1 = pi*(d_fuel**2)/4    #combustion chamber area
if COMBSim:
    T_t = []                #Temperature at the throat
    V_t = []                #Velocity at the throat
    P_t = []                #Pressure at the throat
    vol_1 = []              #Specific volume at the throat
    vol_t = []
    AR = []                 #Area ratio (throat reduction)
    A_t = []                #Throat area
    d_t = []                #Throat diameter
    V_2 = []                #Exhaust velocity
    ISP = []                #Specific impulse
    m_dot = []              #Mass flow
    Thr = []                #Thrust
    burn_rate = []          #Burn rate, in mm/s
    A_comb = []             #Exposed area entering combustion
    d_comb = []             #Diameter of the combustion area (assumed to be the base of a cylinder)
    A_2 = []                #Exit area
    d_2 = []                #Exit diameter       
    for i in range(0,len(P_1)):
        if PressAtm:
            data = combcalc(T_1, k, R, P_1[i], P_atm, desiredThrust, a, n, rho_fuel)
        T_t.append(data[0])
        V_t.append(data[1])
        P_t.append(data[2])
        vol_1.append(data[3])
        vol_t.append(data[4])
        V_2.append(data[5])
        m_dot.append(data[6])
        A_t.append(data[7])
        AR.append(data[8])
        A_2.append(data[9])
        d_t.append(data[10])
        d_2.append(data[11])
        ISP.append(data[12])
        Thr.append(data[13])
        burn_rate.append(data[14])
        A_comb.append(data[15])
        d_comb.append(data[16])              
"""*****************************************************"""
"""NOZZLE CALCULATIONS TABLE"""
if COMBSim:
    colH = ["V_t", "T_t", "P_t", "Pt/P1", "vol_1", "vol_t", "AR", "A_t", "A2", "d_t", "d_2", "V_2", "Isp", "m_dot", "Thrust", "r", "A_c", "d_c", "d_c/d"]
    units = ["m/s", "K", "MPa", "-", "m^3/kg", "m^3/kg", "-", "mm^2", "mm^2", "mm", "mm", "m/s", "s", "kg/s", "N", "mm/s", "mm^2", "mm", "-"]
    data = []
    for i in range(0, n_P):
        data.append(V_t[i])
        data.append(T_t[i])
        data.append(P_t[i]/10**6)
        data.append(P_t[i]/P_1[i])
        data.append(vol_1[i])
        data.append(vol_t[i])
        data.append(AR[i])
        data.append(A_t[i]*10**6)
        data.append(A_2[i]*10**6)
        data.append(d_t[i]*1000)
        data.append(d_2[i]*1000)
        data.append(V_2[i])
        data.append(ISP[i])
        data.append(m_dot[i])
        data.append(Thr[i])
        data.append(burn_rate[i])
        data.append(A_comb[i])
        data.append(d_comb[i])
        data.append(d_comb[i]/d_fuel)
    writetable(fnames[1], colH, P_1, "MPa", 1/(10**6), data, 1, 10, units, "Combustion Dynamics Table")

"""*****************************************************"""

"""*****************************************************"""
"""WALL THICKNESS BASED ON MATERIAL"""
#S. = Seamless (sem costura)
#N1-N2: table number, line number
#Admissible atresses taken at 150 eg C, based on previous results
materials = [["AISI 1020", "AISI 1045", "Glass Fiber", "Alclad 2014-T6", "Carbon Fiber"],   #Name
             [350, 530, 2415, 203, 3000],               #Tensile strength [MPa]
             [7870, 7850, 2110, 2800, 2000],            #Density [kg/m^3]
             [56.7, 56.7, 0.36, 186, 0.628],            #Thermal conductivity [W/mK, 400K]
             [486, 486, 787, 875, 750]]                  #Specific heat [J/kgK]
casingThickness = []
d_inner_casing = []
linCasingMass = []

if not VARMaterials:        #takes only one given material
    temp = materials.copy()
    materials = [[temp[i][chosenMaterial]] for i in range(0,len(temp))]
    temp.clear()

n_Mat = len(materials[0])
for i in range(0, n_P):
    for j in range(0, n_Mat):
        casingThickness.append(casingTCalc(P_1[i], FS, d_casing, materials[1][j]))  #gets casing thickness for every pressure and material
        linCasingMass.append(linMass(d_casing, d_casing - 2*casingThickness[i*n_Mat+j], materials[2][j]))       #gets linear mass of casing for same variables
        d_inner_casing.append(d_casing - 2*casingThickness[i*n_Mat + j])

units = ["mm"]*n_Mat        #thickness in mm
writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), casingThickness, 1000, 16, units, "Casing Thickness")

units = ["kg/m"]*n_Mat        #linear mass in kg/m
writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), linCasingMass, 1, 16, units, "Linear Casing Mass")
"""*****************************************************"""
"""*****************************************************"""
"""TOTAL MASS AND EMPTY MASS CALCULATION"""
if CALCMass:    
    linFuelMass = rho_fuel*pi/4*(d_fuel**2) #fuel occupies the entirety of its diameter
    fuelLength = []
    totalLength = []
    emptyMass = []
    fullMass = []
    for i in range(0, n_P):    #L = -R(mp+mn)/(R(mc + ms) - ms)        
        for j in range(0, n_Mat):
            fuelLength.append(maxHeight)
            # if linFuelMass <= (desiredFuelRatio*linCasingMass[i*n_Mat + j])/(1-desiredFuelRatio):  #this is the situation where no matter how long the fuel grain is, it will never reach the deired fuel/total weight ratio.
            #     fuelLength.append(maxHeight)
            #     print("Configuração {0} MPA e ".format(P_1[i]/10**6) + materials[0][j] + ": fuel ratio inalcançável")
            # else:
            #     tempLength = (desiredFuelRatio*(m_nozzle+m_payload))/(linFuelMass-desiredFuelRatio*(linCasingMass[i*n_Mat + j]+linFuelMass))
            #     if tempLength > maxHeight or tempLength <0:
            #         fuelLength.append(maxHeight)
            #     else:
            #         fuelLength.append(tempLength)
            totalLength.append(l_end +l_payload[2] + fuelLength[i*n_Mat + j])
            emptyMass.append(m_payload + m_nozzle + fuelLength[i*n_Mat + j]*linCasingMass[i*n_Mat + j])
            fullMass.append(emptyMass[i*n_Mat + j] + fuelLength[i*n_Mat + j]*linFuelMass)
    
    fuelMass = [fullMass[m]-emptyMass[m] for m in range(0, len(emptyMass))]
    massRatio = [1-emptyMass[m]/fullMass[m] for m in range(0, len(emptyMass))]  
    
    writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), fuelLength, 1000, 16, "mm"*n_Mat, "Fuel length")
    writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), totalLength, 1000, 16, "mm"*n_Mat, "Total length")

    writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), emptyMass, 1000, 16, ["g"]*n_Mat, "Empty mass")
    writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), fullMass, 1000, 16, ["g"]*n_Mat, "Full mass")
    
    writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), fuelMass, 1000, 16, ["g"]*n_Mat, "Total fuel mass")
    writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), massRatio, 1, 16,  title="Fuel-total mass ratio")


"""*****************************************************"""

"""*****************************************************"""
"""THERMAL CALCULATION"""
if TEMPSim:
    R_cond = []
    R_cond_temp = []
    for i in range(0, n_P):
        for j in range(0, n_Mat):
            R_cond_temp.append(math.log(d_inner_casing[i*n_Mat + j]/d_fuel)/(2*pi*k_paper))
        R_cond.append(R_cond_temp[:])
        R_cond_temp.clear()
"""*****************************************************"""
"""ROCKET TRAJECTORY SIMULATION"""
tlabel = "$t\, [s]$"
xlabel = "$x\, [m]$"
vlabel = "$v\, [m/s]$"
alabel = "$a\, [m/s^2]$"
mlabel = "m\, [kg]"
gllabel = "$G_L\, [m/s]$"
dllabel = "$D_L\, [m/s]$"
tllabel = "$T_L\, [m/s]$"
relabel = "$Re\, [-]$"
malabel = "$Ma\, [-]$"
tslabel ="$T_w\, [K]$"



labels = [tlabel, xlabel, vlabel, alabel, mlabel, gllabel, dllabel, tllabel, relabel, malabel, tslabel]
colors = ['r', 'g', 'b', 'purple']


A_transv = pi/4*d_casing**2
tstep = 0.01                #timestep definition

TWRlist  = []
maxHlist = []

for i in range(0, n_P):
    if CALCMass:
        for j in range(0,n_Mat):
            dotm = 3*fullMass[i*n_Mat+j]/(ISP[i])
            data = trajSim(fullMass[i*n_Mat + j], emptyMass[i*n_Mat + j], ISP[i], dotm, 0.001, d_casing, d_inner_casing[i*n_Mat + j], d_inner_casing[i*n_Mat + j] - 2*t_insulation, materials[2][j], materials[4][j], k_paper,materials[3][j], burn_rate[i])[:]
            maxHlist.append(max(data[1]))
            TWRlist.append(data[3][0]/g0)
            # if TEMPSim:
            #     q = data[-1]                        #gets only the total heat transfer
            #     qt = [sum(q[qi][t] for qi in range(0, len(q))) for t in range(0, len(q[0]))]
            #     data[-1] = qt
            
            for k in range (0, len(data)-1):      #plots all data
                if graphs[k]:                
                    plt.figure(k, dpi=500)
                    plt.xlabel(labels[0])
                    plt.ylabel(labels[k+1])
                    #plt.title("P = {} MPa, FS ={}".format(P_1[i]/1000000, FS))
                    plt.plot(data[0], data[k+1], color=((n_Mat-j)/n_Mat,0, j/n_Mat), linestyle=linstl(j), marker=mrkr(j), markevery = int(10/tstep), label=materials[0][j] +" {:.2f} mm".format(1000*casingThickness[i*n_Mat + j]) + " {:.1f} m".format(totalLength[i*n_Mat + j]))
                    #if i == 0:
                    plt.legend(prop={'size': 7})   
        plt.show()


writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), TWRlist, 1, 16, title="Starting TWR")
writetable(fnames[0], materials[0], P_1, "MPa", 1/(10**6), maxHlist, 1, 16, ["m"]*len(materials[0]), title="Maximum height")
