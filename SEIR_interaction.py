import numpy as np
import matplotlib.pyplot as plt 
from ODESolver import *

class Region:
    def __init__(self, name, S_0, E2_0):
        self.name = name
        self.S_0  = S_0; self.E1_0 = 0; self.E2_0 = E2_0; 
        self.I_0 = 0; self.Ia_0 = 0; self.R_0 = 0;
        
        self.population = sum([self.S_0, self.E2_0,
        self.I_0, self.Ia_0, self.R_0])
       
    def set_SEIR_values(self, u, t):
       self.S = u[:,0]; self.E1 = u[:,1]; self.E2 = u[:,2];
       self.I = u[:,3]; self.Ia = u[:,4]; self.R = u[:,5]
       self.t = t
    
    
    def plot(self):
        S = self.S; I = self.I; 
        Ia = self.Ia; R = self.R;
        t = self.t
        
        plt.plot(t, S, label = "Susceptible")
        plt.plot(t, I, label = "Symtomatically infected ")
        plt.plot(t, Ia, label = "Asymptomatically infected")
        plt.plot(t, R, label = "Resistant")
        plt.title(self.name)
        plt.xlabel("Time(days")
        plt.ylabel("Population")
        
class RegionInteraction(Region):
    
    def __init__(self, name, S_0, E2_0, lat, long):
        self.S_0  = S_0; self.E1_0 = 0; self.E2_0 = E2_0; 
        self.I_0 = 0; self.Ia_0 = 0; self.R_0 = 0;
        
        self.lat = lat*(np.pi/180)
        self.long = long*(np.pi/180)
        self.name = name
        
        self.population = sum([self.S_0, self.E2_0,
        self.I_0, self.Ia_0, self.R_0])
        
        
    def distance(self, other):
        
        lat = self.lat; long = self.long
        R_earth = 64
        phi_i = lat; phi_j = other.lat;
        lambda_i = long; lambda_j = other.long
        arg = np.sin(phi_i)*np.sin(phi_j)+np.cos(phi_i)*np.cos(phi_j)*np.cos(abs(lambda_i-lambda_j))
        if not 0<=abs(arg)<=1:
            print(f"the argument: {arg} is not between 0 and 1")
        sigma_ij = np.arccos(arg)
        d = R_earth*sigma_ij
        return d 

class ProblemSEIR():
    def __init__(self, region, beta, r_ia = 0.1, r_e2=1.25,lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2):
        if isinstance(beta, (float, int)): 
            self.beta = lambda t: beta    
        elif callable(beta):
            self.beta = beta
        self.region = region
        self.r_ia = r_ia; self.r_e2 = r_e2; 
        self.lmbda_1 = lmbda_1; self.lmbda_2 = lmbda_2;
        self.p_a = p_a; self.mu = mu
        self.set_initial_condition()
        
    def set_initial_condition(self):
        S_0 = self.region.S_0
        E2_0 = self.region.E2_0
        self.initial_condition = [S_0, 0, 
                                  E2_0,
                                  0, 0, 0]
        
    def get_population(self):
        return self.region.population
            
    def solution(self, u, t):
        self.region.set_SEIR_values(u,t)
        
        
    def __call__(self, u, t):
         beta = self.beta(t); r_ia = self.r_ia; r_e2=self.r_e2;
         lmbda_1=self.lmbda_1; lmbda_2=self.lmbda_2;
         p_a=self.p_a; mu=self.mu;
         S, E1, E2, I, Ia, R = u
         N =sum(u)
         dS  = -beta*S*I/N - r_ia*beta*S*Ia/N - r_e2*beta*S*E2/N
    
         dE1 = beta*S*I/N + r_ia*beta*S*Ia/N + r_e2*beta*S*E2/N - lmbda_1*E1
         dE2 = lmbda_1*(1-p_a)*E1 - lmbda_2*E2
         dI  = lmbda_2*E2 - mu*I
         dIa = lmbda_1*p_a*E1 - mu*Ia
         dR  = mu*(I + Ia)
         return[dS, dE1, dE2, dI, dIa, dR]
     
class ProblemInteraction(ProblemSEIR):
    def __init__(self, region, area_name, beta, r_ia = 0.1, r_e2=1.25,lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2):
        if isinstance(beta, (float, int)): 
            self.beta = lambda t: beta    
        elif callable(beta):
            self.beta = beta
        self.area_name = area_name
        self.region = region
        self.r_ia = r_ia; self.r_e2 = r_e2; 
        self.lmbda_1 = lmbda_1; self.lmbda_2 = lmbda_2;
        self.p_a = p_a; self.mu = mu
        self.set_initial_condition()
        
    def get_population(self):
        region = self.region
        pop = []
        for k in range(len(region)):
            val = region[k].population
            pop.append(val)
        return sum(pop)
    
    def set_initial_condition(self):
        region = self.region
        M = len(region)
        self.initial_condition = []
        for i in range(M):
            self.initial_condition += [region[i].S_0, 
                                      region[i].E1_0, 
                                      region[i].E2_0, 
                                      region[i].I_0,
                                      region[i].Ia_0,
                                      region[i].R_0]
        
    def __call__(self, u, t):
        beta = self.beta(t); r_ia = self.r_ia; r_e2=self.r_e2;
        lmbda_1=self.lmbda_1; lmbda_2=self.lmbda_2;
        p_a=self.p_a; mu=self.mu;
         
        n = len(self.region)
        
        SEIR_list = [u[i:i+6] for i in range(0, len(u), 6)]
        E2_list = [u[i] for i in range(2, len(u), 6)]
        Ia_list = [u[i] for i in range(4, len(u), 6)]
        derivative = []
        for i in range(n):
            S, E1, E2, I, Ia, R = SEIR_list[i]
            dS = 0; 
            N_i = self.region[i].population
            
            for j in range(n):
                E2_other = E2_list[j]
                Ia_other = Ia_list[j]
                N_j = self.region[j].population
                dij = self.region[i].distance(self.region[j])
               
                dS  += (-beta*S*I)/N_i - (r_ia*beta*S)*((Ia_other/N_j)*np.exp(-dij)) - r_e2*beta*S*((E2_other/N_j)*np.exp(-dij))
           
            dE1 = -dS - lmbda_1*E1
            dE2 = lmbda_1*(1-p_a)*E1 - lmbda_2*E2
            dI  = lmbda_2*E2 - mu*I
            dIa = lmbda_1*p_a*E1 - mu*Ia
            dR  = mu*(I + Ia)
            
            derivative += [dS, dE1, dE2, dI, dIa, dR]
        
        return derivative
                
                
    def solution(self, u, t):
            n = len(t)
            n_reg = len(self.region)
            self.t = t
            self.S = np.zeros(n)
            self.E1 = np.zeros(n)
            self.E2 = np.zeros(n)
            self.I = np.zeros(n)
            self.Ia = np.zeros(n)
            self.R = np.zeros(n)
            SEIR_list = [u[:, i:i+6] for i in range(0, n_reg*6, 6)]
            for part, SEIR in zip(self.region, SEIR_list):
                part.set_SEIR_values(SEIR, t)
                self.S += SEIR[:,0]
                self.E1 += SEIR[:,1]
                self.E2 += SEIR[:,2]
                self.I += SEIR[:,3]
                self.Ia += SEIR[:,4]
                self.R += SEIR[:,5]
        
    
    def plot(self):
        S = self.S; I = self.I; 
        Ia = self.Ia; R = self.R;
        t = self.t
       
        plt.plot(t, S, label = "Susceptible")
        plt.plot(t, I, label = "Symtomatically infected ")
        plt.plot(t, Ia, label = "Asymptomatically infected")
        plt.plot(t, R, label = "Resistant")
        plt.title(self.area_name)
        plt.xlabel("Time(days")
        plt.ylabel("Population")
        
         
            
            
        
        
class SolverSEIR:
    def __init__(self, problem, T, dt):
        self.problem = problem
        self.T = T
        self.dt = dt
        self.total_populatiion = problem.get_population
        
    def solve(self, method = RungeKutta4):
        solver = method(self.problem)
        solver.set_initial_condition(problem.initial_condition)
        t = np.linspace(0,self.T, self.T*int(self.dt))
        u, t = solver.solve(t)
        self.problem.solution(u, t)
        
        
        
if __name__== "__main__":
    innlandet = RegionInteraction("Innlandet",S_0=371385, E2_0=0, lat=60.7945,long=11.0680)
    oslo = RegionInteraction("Oslo",S_0=693494,E2_0=100, lat=59.9,long=10.8)
    print(oslo.distance(innlandet))
    print("-----------------------------------")
    problem = ProblemInteraction([oslo, innlandet],"Norway_east", beta=0.5)
    print(problem.get_population())
    problem.set_initial_condition()
    print(problem.initial_condition)
    u = problem.initial_condition

    print(problem(u,0))
    solver = SolverSEIR(problem,T=100,dt=1.0)
    solver.solve()
    problem.plot()
    plt.legend()
    plt.show()

"""
% python3 SEIR_interaction.py
1.0100809386285283
-----------------------------------
1064979
[693494, 0, 100, 0, 0, 0, 371385, 0, 0, 0, 0, 0]
[-62.490904683238746, 62.490904683238746, -50.0, 50.0, 0.0, 0.0, -12.1878323242723, 12.1878323242723, 0.0, 0.0, 0.0, 0.0]

"""