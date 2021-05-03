import numpy as np
import matplotlib.pyplot as plt
import random

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

def get_C0(gamma,Beta):
    """Given the current gamma and beta, returns the updated C0 value
    C0 is defined as the right side of the Bethe formula divided by rho"""
    C0=0.3071*Z*z*z/(A*Beta*Beta)*(np.log(2*m*c*c*Beta*Beta*gamma*gamma/I)-Beta*Beta)
    return C0

def get_rho(x):
    """Given the height from the ground, returns the mass dencity, rho of the air at that height."""
    rho=p0*M/(R*T0)*(1-L*x/T0)**(g*M/(R*L)-1)
    return rho

#Hey group members can yall check which of these get_t_prime functions is right????? - thx

# def alternate_get_t_prime(C0, rho, gamma1, gamma2):
#     """Given the updated C0 and rho, the previous gamma, and the updated gamma, returns the updated t_prime
#     t_prime is the time elapsed in the muon's reference frame"""
#     t_prime = m*c/(rho*C0)*(gamma1 * (gamma1**2 - 1)**(-3/2) - gamma2 * (gamma2**2 - 1)**(-3/2))
#     return t_prime

# def get_t_prime(C0,rho,gamma1,gamma2):
#     t_prime=m*c/(rho*C0)*(1/np.sqrt(gamma2**2-1)-1/np.sqrt(gamma1**2-1))*(gamma2-gamma1)
#     return t_prime


# we're just doing riemann sum here, I think. So it should be the value of the expression, times dgamma
def get_dt_prime(C0,rho,gamma1,gamma2):
	gamma_average = (gamma2+gamma1)/2.
	dg = gamma1-gamma2
	print(C0)
	dt_prime = (m*c/(rho*C0)) * dg * 1./(np.sqrt(gamma_average**2-1))
	return dt_prime

def findDecayProbability(t_prime):
    """Given the the time it takes for the muon to reach the detector in its own reference frame, returns the probability that it decays. """
    tau = 2.2*10**-6 #seconds
    return np.exp(-1*t_prime/tau)
    

#Initial Conditions, will be provided from the previous groups as a list for each condition of all the muons
E_initial = np.array([1e18,1e20], dtype = 'f')
x0_flat_initial = np.array([12000,10000], dtype = 'f') #height of troposphere
x0_round_initial = np.array([15000,13000], dtype = 'f')

#instead of appending, it makes the list of intial gammas by using numpy. A litte faster - mason
gamma_initial = E_initial/(m*c**2)

#same thing here
Beta_initial = np.sqrt(1-1/(gamma_initial**2))

t_prime_flat = 0
t_prime_round = 0


E_flat_final = np.zeros(len(x0_flat_initial))
E_round_final = np.zeros(len(x0_flat_initial))
t_prime_flat_final = np.zeros(len(x0_round_initial))
t_prime_round_final = np.zeros(len(x0_round_initial))

number_of_steps = 1000. #dummy amount, will change
dx_flat = x0_flat_initial/number_of_steps
dx_round = x0_round_initial/number_of_steps

#Find final energies and t_prime for flat earth
Beta = Beta_initial
x = x0_flat_initial
gamma = gamma_initial
E = E_initial

for i in range(int(number_of_steps)):
    C0=get_C0(gamma,Beta)
    rho=get_rho(x)
    dE=C0*rho*dx_flat
    E1 = E
    E2 = E1 - dE #feel like this should be adding dE but that results in increasing energy
    gamma1 = gamma
    gamma2 = E2/(m*c**2)
    Beta1=Beta
    Beta2=np.sqrt(1-1/(gamma2*gamma2))
    t_prime_flat += get_dt_prime(C0,rho,gamma1,gamma2)

    gamma = gamma2
    Beta = Beta2
    E = E2
    x += dx_flat

E_flat_final = E
t_prime_flat_final = t_prime_flat


#Find the final energies and t_prime for round earth
Beta = Beta_initial
x = x0_round_initial
gamma = gamma_initial
E = E_initial

for i in range(int(number_of_steps)):
    C0=get_C0(gamma,Beta)
    rho=get_rho(x)
    dE=C0*rho*dx_round
    E1 = E
    E2 = E1 - dE #feel like this should be adding dE but that results in increasing energy
    gamma1 = gamma
    gamma2 = E2/(m*c**2)
    Beta1=Beta
    Beta2=np.sqrt(1-1/(gamma2*gamma2))
    t_prime_round += get_dt_prime(C0,rho,gamma1,gamma2)

    gamma = gamma2
    Beta = Beta2
    E = E2
    x += dx_round

E_round_final = E
t_prime_round_final = t_prime_round
 

print(t_prime_round)
print(t_prime_flat)

#Simulating how many muons hit the detector:
    #Checks for muons in the proper energy range
    #Makes a probablistic decision whether a paticular muon has decayed
    #If muon has the proper energy and is found not to decay, add to histogram
        
energyDetected = 160 #Mev
energyAllowance = 2 #Mev
 
anglesOfDetectedMuons_flatEarth = []
anglesOfDetectedMuons_roundEarth = []

for muon in range(len(x0_flat_initial)):
    #Check detection for flat Earth
    if (energyDetected - energyAllowance <= E_flat_final[muon] <= energyDetected + energyAllowance):
        #Decide if it decayed
        if (random.random() < findDecayProbability(t_prime_flat_final[muon])):
            anglesOfDetectedMuons_flatEarth.append(zenithAngle_final[muon])
        
    #Check detection for round Earth
    if (energyDetected - energyAllowance <= E_round_final[muon] <= energyDetected + energyAllowance):
        #Decide if it decayed
        if (random.random() < findDecayProbability(t_prime_round_final[muon])):
            anglesOfDetectedMuons_roundEarth.append(zenithAngle_final[muon])
        
bin_number = 20       

plt.hist(anglesOfDetectedMuons_flatEarth, bins = bin_number)
plt.title('Flat Earth Distribution of Detected Muons as a Function of Zenith Angle')
plt.xlabel('Zenith Angle (Radians)')
plt.ylabel('Number of Detected Muons')
plt.show()

plt.hist(anglesOfDetectedMuons_roundEarth, bins = bin_number)
plt.title('Round Earth Distribution of Detected Muons as a Function of Zenith Angle')
plt.xlabel('Zenith Angle (Radians)')
plt.ylabel('Number of Detected Muons')
plt.show()
