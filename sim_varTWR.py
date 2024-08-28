import matplotlib.pyplot as plt
import math

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

g0 = 9.81   #CONSTANT
Isp = 130   #assumed constant (for now); will only ever increase during ascent

m0 = 3.0    #rocket mass parameters
me = 1.0

Cd = 0.5    #entirely guessed, needs further research
d = 0.1
A = 3.14*(d**2)/4       #supposing 10 cm diameter rocket

TWR_low = 1.0
TWR_high = 5.0
TWR_step = 0.5
n_TWR = int((TWR_high-TWR_low)/TWR_step)+1


tstep = 0.05         #timestep definition
        
h = []      #atmospheric properties
#Tp = []
g = []
P = []
rho = []
mu = []
c = 330                     #assumed constant (about right until 20km)

def Re(rho, v, d, mu):      #calculates Reynolds' number
    return (rho*v*d/mu)

def Ma(v):                  #assumes constant speed of sound (fair)
    return (v/c)

def linintrp (xaxis, yaxis, xinterp):       #linear interpolation
    for x in range(0, len(yaxis)):
        if xinterp <=0:
            return yaxis[0]
        if (xinterp >= xaxis[x] and xinterp <xaxis[x+1]):
            x1 = xaxis[x]
            y1 = yaxis[x]
            x2 = xaxis[x+1]
            y2 = yaxis[x+1] 
            return (y1 + (y2-y1)/(x2-x1)*(xinterp-x1))

f = open("properties.txt", "r")
n_props = len(f.readlines()) - 1
f.close()
f = open("properties.txt", "r")
f.readline()
f.readline()
for line in f:
    spl = line.split("\t")
    h.append(float(spl[0]))
    #Tp.append(float(spl[1]))
    g.append(float(spl[2]))
    P.append(float(spl[3])*10**4)
    rho.append(float(spl[4]))
    mu.append(float(spl[5])*10**-5)
f.close()
max_h = 0.0
max_h_TWR = 0

for i in range (0, n_TWR):
    
    TWR = [TWR_low + TWR_step*i]         #initial thrust-to-weight-ratio
    m = [3.0]           #initial mass
    W = [g0*m[0]]       #initial weight
    T = [TWR[0]*W[0]]   #initial thrust
    D = [0.0]           #initial drag
    FR = [T[0]-W[0]-D[0]]
    Re_t = [0.0]
    Ma_t = [0.0]
    mdot = T[0]/(Isp*g0)    #mass consumption rate
    
    deltaV = [Isp*g0*math.log(m[0])/me]
    
    time = [0.0]        #initial time
    
    a = [FR[0]/m[0]]       #initial acceleration
    v = [0.0]           #initial speed
    x = [0.0]           #initial position
    g_loss = [0.0]
    D_loss = [0.0]
    tot_loss = [0.0]
    t_inv = 0.0
    
    t=0
    while x[t]>=0.0:
        time.append(time[t]+tstep)
        if (m[t]<=me):  #i.e. if the rocket ran out of fuel
            T.append(0)
            m.append(m[t])
            g_loss.append(g_loss[t])
        else:
            T.append(T[0])           #if i ever do change the thrust calculation, this may be handy
            m.append(m[t]- (mdot*tstep))
            g_loss.append(g_loss[t]+(W[t]/T[t])*g0*tstep)
        
        W.append(m[t]*linintrp(h, g, x[t]))
        
        if v[t]>=0:     #if going up
            D.append(linintrp(h, rho, x[t])*A*Cd*v[t]**2/2)     #drag is in -x
            D_loss.append(D_loss[t]+D[t]*tstep/m[t])
        else:
            D.append(-linintrp(h, rho, x[t])*A*Cd*v[t]**2/2)    #drag is in +x
            D_loss.append(D_loss[t])
        FR.append(T[t]-W[t]-D[t])
        a.append(FR[t]/m[t])
        v.append(v[t]+a[t]*tstep)
        x.append(x[t]+v[t]*tstep)
            
        TWR.append(T[t]/W[t])
        Re_t.append(Re(linintrp(h, rho, x[t]), v[t], d, linintrp(h, mu, x[t])))
        Ma_t.append(Ma(v[t]))
        
        tot_loss.append(g_loss[t]+D_loss[t])
        deltaV.append(Isp*g0*math.log(m[t])/me)
        if D > W:
            t_inv = t*tstep
        t +=1
        
    for j in range(0, len(x)):      #gets maximum height and twr to get it
        if x[j] > max_h:            
            max_h = x[j]
            max_h_TWR = i
    
    TWRlabel = str(TWR[0])
    TWRlabel = TWRlabel.replace('.', ',')        
    
    plt.figure(0, dpi=500)
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$x\,[m]$")
    plt.plot(time, x, color=((n_TWR-i)/n_TWR,0,i/n_TWR), marker=mrkr(i), markevery=int(10/tstep), linestyle=linstl(i), label="TWR {0}".format((TWRlabel)))
    plt.legend(prop={'size': 8})
    
    plt.figure(1, dpi=500)
    #plt.title("V x t")
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$V\,[m/s]$")
    plt.plot(time, v, color=((n_TWR-i)/n_TWR,0,i/n_TWR), marker=mrkr(i), markevery=int(10/tstep), linestyle=linstl(i), label="TWR {0}".format((TWRlabel)))
    #plt.legend(prop={'size': 8})

    plt.figure(2, dpi=500)
    #plt.title("a x t")
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$a\,[m/s^2]$")
    plt.plot(time, a, color=((n_TWR-i)/n_TWR,0,i/n_TWR), marker=mrkr(i), markevery=int(10/tstep), linestyle=linstl(i), label="TWR {0}".format((TWRlabel)))
    #plt.legend(prop={'size': 8})
    
    plt.figure(3, dpi=500)
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$G_L\,[m/s]$")
    plt.plot(time, g_loss, color=((n_TWR-i)/n_TWR,0,i/n_TWR), marker=mrkr(i), markevery=int(10/tstep), linestyle=linstl(i), label="TWR {0}".format((TWRlabel)))
    #plt.legend(prop={'size': 8})
    
    plt.figure(4, dpi=500)
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$D_L\,[m/s]$")
    plt.plot(time, D_loss, color=((n_TWR-i)/n_TWR,0,i/n_TWR), marker=mrkr(i), markevery=int(10/tstep), linestyle=linstl(i), label="TWR {0}".format((TWRlabel)))
    #plt.legend(prop={'size': 8})
    
    plt.figure(5, dpi=500)
    #plt.title("Total losses x t")
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$T_L\,[m/s]$")
    plt.plot(time, tot_loss, color=((n_TWR-i)/n_TWR,0,i/n_TWR), marker=mrkr(i), markevery=int(10/tstep), linestyle=linstl(i), label="TWR {0}".format((TWRlabel)))
    #plt.legend(prop={'size': 8})
    
    plt.figure(6, dpi=500)
    #plt.title("DeltaV x t")
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$\Delta V\,[m/s]$")
    plt.plot(time, deltaV, color=((n_TWR-i)/n_TWR,0,i/n_TWR), marker=mrkr(i), markevery=int(10/tstep), linestyle=linstl(i), label="TWR {0}".format((TWRlabel)))
    #plt.legend(prop={'size': 8})
    
    #resets all lists    
    time.clear()
    a.clear()
    v.clear()
    x.clear()
    m.clear()
    TWR.clear()
    W.clear()
    T.clear()
    D.clear()
    FR.clear()
    Re_t.clear()
    Ma_t.clear()
    m.clear()
 
plt.show()  
        
print("m0: {0}\tmf: {1}\t m0/mf: {2}".format(m0, me, m0/me))
print("Cd: {0}\td: {1}\t A: {2}".format(Cd, d, 3.14*d**2/4))
print("Isp: {0}".format(Isp))
print("TWR de máxima altura: {0}\nAltura máxima: {1}".format(TWR_low + max_h_TWR*TWR_step, max_h))