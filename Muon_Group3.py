import numpy as np

#Constants: Not subject to change
c=3e8 #Speed of light, m/s
m=1.88e-28/1.6605e-27 #mass of muon, AU
p0=101325 #sea level standard atmospheric pressure, Pa
T0= 288.15 #sea level standard temperature, K
g= 9.807 #m/s^2
L= .0065 #temperature lapse rate, K/m
R= 8.31446 #ideal gas constant, J/(mol K)
M= 0.0289652 #molar mass of dry air, kg/mol
e = -1.602e-19 #electron charge, e

#Constants: Subject to change based on info we get from previous groups
A=13 #atomic weight of the stopping medium
Z=7 #atomic number of N
I=A*Z #approximation of mean excitation energy of the stopping medium of atoms
z=1*e #charge of muon relative to electron charge

def get_C0(x,gamma,Beta):
    """Given the current height from the ground, gamma, and beta, returns the updated C0 value
    C0 is defined as the right side of the Bethe formula divided by rho"""
    C0=0.3071*Z*z*z/(A*Beta*Beta)*(np.log(2*m*c*c*Beta*Beta*gamma*gamma/I)-Beta*Beta)
    return C0

def get_rho(x):
    """Given the height from the ground, returns the mass dencity, rho of the air at that height."""
    rho=p0*M/(R*T0)*(1-L*x/T0)**(g*M/(R*L)-1)
    return rho

#Hey group members can yall check which of these get_t_prime functions is right????? - thx

def alternate_get_t_prime(C0, rho, gamma1, gamma2):
    """Given the updated C0 and rho, the previous gamma, and the updated gamma, returns the updated t_prime
    t_prime is the time elapsed in the muon's reference frame"""
    t_prime = m*c/(rho*C0)*(gamma1 * (gamma1**2 - 1)**(-3/2) - gamma2 * (gamma2**2 - 1)**(-3/2))
    return t_prime

def get_t_prime(C0,rho,gamma1,gamma2):
    t_prime=m*c/(rho*C0)*(1/np.sqrt(gamma2**2-1)-1/np.sqrt(gamma1**2-1))*(gamma2-gamma1)
    return t_prime
    

#Initial Conditions, will be provided from the previous groups as a list for each condition of all the muons
Beta_initial=[.1,.1] #dummy variable amount, v/c
x0_initial=[12000,10000] #height of troposphere
gamma_initial = [] #Lorenz Factor
for j in range(len(Beta_initial)):
    gamma_initial.append(1/np.sqrt(1-Beta_initial[j]*Beta_initial[j]))
E_initial=[] #double check, not correct
for k in range(len(x0_initial)):
    E_initial.append(gamma_initial[k]*m*c**2)
t_prime = 0

deltax=1 #dummy amount, will change

Beta_final = []
E_final = []
C0_final = []
rho_final = []
gamma_final = []
t_prime_final = []

for muon in range(len(x0_initial)):
    Beta = Beta_initial[muon]
    x0 = x0_initial[muon]
    gamma = gamma_initial[muon]
    E = E_initial[muon]
    for i,x in enumerate(range(x0,0,-deltax)):
        C0=get_C0(x,gamma,Beta)
        rho=get_rho(x)
        dE=C0*rho*deltax
        E1 = E
        E2 = E1 - dE #feel like this should be adding dE but that results in increasing energy
        gamma1 = gamma
        gamma2 = E/(m*c**2)
        #gamma2=E2*gamma1/E1
        Beta1=Beta
        Beta2=np.sqrt(1-1/(gamma2*gamma2))
        t_prime =+ get_t_prime(C0,rho,gamma1,gamma2)

        gamma = gamma2
        Beta = Beta2
        E = E2
        if x%500==0:
            print(x,Beta,E,C0,rho,gamma,t_prime,'\n')
    Beta_final.append(Beta)
    E_final.append(E)
    C0_final.append(C0)
    rho_final.append(rho)
    gamma_final.append(gamma)
    t_prime_final.append(t_prime)



