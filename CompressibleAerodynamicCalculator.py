import numpy as np
import math
from math import sqrt, sin, cos, tan, atan
from scipy.optimize import fsolve

class IsentropicFlowRelations():
    def __init__(self, gamma, R, Mach):
        self.gamma = gamma
        self.R = R
        self.Mach = Mach

    def staticByTotalPressure(self):
        return pow(1 +  (self.gamma - 1)/2 * self.Mach**2, -(self.gamma/(self.gamma-1)))

    def staticByTotalTemperature(self):
        return pow(1 +  (self.gamma - 1)/2 * self.Mach**2, -1)

    def staticByTotalDensity(self):
        return pow(1 +  (self.gamma - 1)/2 * self.Mach**2, (-1/(self.gamma-1)))

    def staticByChokeArea(self):
        return  pow(1 + self.Mach**2 * (self.gamma - 1)/2 , (self.gamma + 1)/(2*(self.gamma-1))) * pow((self.gamma + 1)/2, -((self.gamma + 1)/(2*(self.gamma -1))))/self.Mach

    def machAngle(self): #degree
        print("Mach Angle is in Degree.")
        return math.asin(1/self.Mach)* (180/math.pi)

    def prandtlMeyerAngle(self):
        a = sqrt((self.gamma + 1) / (self.gamma - 1)) * math.atan(sqrt((self.gamma - 1) / (self.gamma + 1) * (self.Mach**2 - 1)))
        b = math.atan(sqrt(self.Mach**2 - 1))
        return (a - b) * (180 / math.pi)

class NormalShockRelations(IsentropicFlowRelations):
    def secondMach(self):
        a = ((self.gamma-1)* self.Mach**2) + 2
        b = 2*self.gamma* self.Mach**2 - (self.gamma-1)
        return sqrt(a/b)

    def AfterDensity(self):
        return ((self.gamma + 1)* self.Mach**2)/ ((self.gamma-1)* self.Mach**2 + 2)

    def AfterTemperature(self):
        a = (2* self.Mach**2* self.gamma - (self.gamma -1)) * ((self.gamma-1)* self.Mach**2 + 2)
        b = (self.gamma + 1)**2 * self.Mach**2
        return a/b
    def AfterPressure(self):
        a = 2*self.gamma*self.Mach**2 - (self.gamma - 1)
        b = self.gamma + 1
        return a/b

    def AfterTotalPressure(self):
        a = pow(self.AfterDensity(), self.gamma/(self.gamma -1))
        b = pow(1/self.AfterPressure(), 1/(self.gamma -1))
        return a*b

class ObliqueShockRelations(IsentropicFlowRelations):
    def __init__(self,gamma, R, Mach, turnAngle):
        super().__init__(gamma, R, Mach)
        # converting the degree into Radians
        self.turnAngle = np.deg2rad(turnAngle)

    def theta_beta_m_relation(self, beta):
        """Theta-Beta-Mach relation equation."""
        beta = np.deg2rad(beta)
        term1 = 2 * (1 / np.tan(beta)) * (self.Mach**2 * np.sin(beta)**2 - 1)
        term2 = self.Mach**2 * (self.gamma + np.cos(2 * beta)) + 2
        return np.tan(self.turnAngle) - (term1 / term2)

    #TODO Find shock angle (beta) using numerical solver.
    def find_shock_angle(self):
        # Initial guess in degree
        beta_guess = 45
        beta_solution = fsolve(self.theta_beta_m_relation, beta_guess)
        return beta_solution[0]

    def deflection_angle(self):
        beta = self.find_shock_angle()
        return np.rad2deg(self.turnAngle)

    def changeInDensity(self):
        a = (self.gamma + 1)* self.Mach**2* sin(np.deg2rad(self.find_shock_angle()))**2
        b = (self.gamma - 1)* self.Mach**2* sin(np.deg2rad(self.find_shock_angle()))**2 + 2
        return (a/b)

    def changeInPressure(self):
        a = 2* self.gamma * self.Mach**2 * sin(np.deg2rad(self.find_shock_angle()))**2 - (self.gamma - 1)
        return a/(self.gamma + 1)

    def changeInTemperature(self):
        return self.changeInPressure() * (1/self.changeInDensity())

    def totalPressureChange(self):
        a = pow(self.changeInDensity(), self.gamma/(self.gamma - 1))
        b = pow(1/self.changeInPressure(), 1/(self.gamma - 1))
        return a*b

    def AfterShockMachNumber(self):
        a = 2 + (self.gamma -1) * self.Mach**2 * (sin(np.deg2rad(self.find_shock_angle())))**2
        b = 2 * self.gamma * self.Mach**2 * (sin(np.deg2rad(self.find_shock_angle())))**2 - (self.gamma -1)
        c = np.sqrt(a/b)
        return (1/(sin(np.deg2rad(self.find_shock_angle()) - (self.turnAngle)))) * (c)


p1 = IsentropicFlowRelations(1.4, 343.24, 2)
p2 = NormalShockRelations(1.4, 343.24, 2)
p3 = ObliqueShockRelations(1.4, 343.24, 2, 20)
print("P/Po: %f" % p1.staticByTotalPressure())
print("T/To: {}".format(p1.staticByTotalTemperature()))
print("rho/rho_o: %f" % p1.staticByTotalDensity())
print("A/A*: %f" % p1.staticByChokeArea())
print("Mach Angle: %f" % p1.machAngle())
print("P-M Angle: %f" % p1.prandtlMeyerAngle())
print("After Normal Shock M: %f" % p2.secondMach())
print("Density change: %f" % p2.AfterDensity())
print("Temperature change: %f" % p2.AfterTemperature())
print("Pressure change: %f" % p2.AfterPressure())
print("Total Pressure change: %f" % p2.AfterTotalPressure())
print("Shock Angle: %f" % p3.find_shock_angle())
print("Change in Density: %f" % p3.changeInDensity())
print("Change in Pressure: %f" %p3.changeInPressure())
print("Change in Temperature: %f" % p3.changeInTemperature())
print("Total Pressure Change: %f" % p3.totalPressureChange())
print("After Shock Mach Number: %f" % p3.AfterShockMachNumber())


