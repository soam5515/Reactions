import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import special
from itertools import product
from scipy.optimize import fsolve, curve_fit
from scipy.interpolate import UnivariateSpline

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
def findpeakinfo(xs,arr):
    spline = UnivariateSpline(xs, [a-np.max(arr)/2 for a in (arr)], s=0)
    print(spline)
    print(spline.roots())
    r0, r1, r2 = spline.roots() # fi    nd the roots
    return np.max(arr), np.abs(r1-r2)

def makeanglecontinuous(angles):
    "Renormalize Angles to remove discontinuities by adding n*pi"
    count=0#how many times n pi we need to add/subtract
    outarray=[]
    for theta in angles:
        if(theta>=0):
            outarray.append(theta)
        else:
            outarray.append(theta+np.pi)
    return outarray

def diffeq(L,Energy):
    "This function returns a function that will be the differential equation to be solved based on E and L"
    def deriv_chi(f, r):
            chi, chiprime = f
            return [chiprime, (-2.9233295/(1.0+np.exp((r-2.58532)/.65))*chi-.047845*Energy*chi+L*(L+1)*chi/r**2) ]
    return deriv_chi
sin=np.sin
def HenkelPlus(rho,L):#Defining Hankel functions and their derivatives
    return np.exp(1j*(rho-L*np.pi/2))
def HenkelMinus(rho,L):
    return np.exp(-1j*(rho-L*np.pi/2))
def HenkelPlusPrime(rho,L,k):
    return k*1j*np.exp(1j*(rho-L*np.pi/2))
def HenkelMinusPrime(rho,L,k):
    return k*(-1j)*np.exp(-1j*(rho-L*np.pi/2))



class chi:
    "At some level a general second order differential equation solver, which anticipates the problem. It uses runge-kutta in SolveDifEq, after using SetDiffEq to make the equation to be solved. It keeps the result, as well as all input used in an object for easy calculation and plotting"
    diff=0#Differential equation to be solved, set via SetDiffEq
    init=[]#array of initial values in the form [y(0),y'(0)]
    rvalues=[]#array containing all r values to be sampled over 
    Energy=0#Energy of system in Schrodinger Eq.
    L=0#Angular Momentum Quantum number in Schro. Eq.
    chivalues=[]#Array containing chi(x) values, each element in the array corresponds as such: chivalues[i]=chi(rvalues[i])
    chiprimevalues=[]#Similar to chivalues but for chi prime
    delta=0#Phase shift.
    Rmat=0#Logarithmic Derivative divided by a
    SMat=0#Represents admixture of Hankel plus function(which is to say "outgoing spherical wave")
    sindelta=0#sin of delta
    Mass=931.49272#Where this comes from taking .0478 from the notes, dividing by 2, multipling by hbar, and then dividing by (1 fm^2*1MeV)
    crosssection=0#Cross Section in arb units.
    aindex=0#index of chivalues, chiprimevalues to be used to generate R matrix(i.e. chivalues[aindex]/chivaluesprime[aindex])
    def SetDiffEq(self):
        "Choose Differential Equation to be solved"
        self.diff=diffeq(self.L,self.Energy)
    def SolveDifeq(self):
        "Solve Differential Equation"
        z = integrate.odeint(self.diff, self.init, self.rvalues,full_output=0,mxstep=10000)
        self.chivalues, self.chiprimevalues=z.T
    def GetK(self):
        "Calculate k in fm"
        k=np.sqrt(.047845*self.Energy)#hbar^2k^2/2m=E
        return k
    def SetRMatrix(self):
        "Find Logarithmic Derivative"
        self.Rmat=self.chivalues[self.aindex]/self.chiprimevalues[self.aindex]/self.rvalues[self.aindex]
    def SetSMatrix(self):
        "Calculate S Matrix from RMatrix using Hankel Functions"
        k=self.GetK()
        a=self.rvalues[self.aindex]
        num=HenkelMinus(k*a,self.L)-a*self.Rmat*HenkelMinusPrime(k*a,self.L,k)
        denom=HenkelPlus(k*a,self.L)-a*self.Rmat*HenkelPlusPrime(k*a,self.L,k)
        self.SMat=num/denom
    def SetDelta(self):
        "Calculate Delta From Logarithm"
        self.delta=1/(2j)*np.log(self.SMat)
    def SetSinDelta(self):
        "Calculate Sin Delta, save in chi"
        self.sindelta=(sin(np.real(self.delta)))
    def SetCrossSection(self):
        "Calculate Cross Section in units of barns"
        self.crosssection=4*np.pi*(2*self.L+1)/(100*self.GetK()**2)*(self.sindelta)**2
    def __init__(self, rvalues,Energy,L,init,aindex):
        "Initialization, requires range, Energy, L, initial conditions and the index to be sampled at."
        self.rvalues=rvalues
        self.Energy=Energy
        self.L=L
        self.init=init
        self.aindex=aindex
def twodplot(x,y,title,xaxis,yaxis):
    fig=plt.figure()
    ax2=fig.add_subplot(111)
    ax2.plot(x,y)
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
        

init=[.001,.001]
r=np.linspace(.0001,100,1000)
Energies=[.1,10]   
Ls=[0,1,2]
radwaves=[]
ain=-10
for vars in product(Energies,Ls):
    radwav=chi(r,vars[0],vars[1],init,ain)
    radwav.SetDiffEq()
    radwav.SolveDifeq()
    radwav.SetRMatrix()
    radwav.SetSMatrix()
    radwav.SetDelta()
    radwav.SetSinDelta()
    twodplot(radwav.rvalues,radwav.chivalues,"$\chi$ vs r at Energy "+str(vars[0])+" and L "+str(vars[1]),"r (fm)", r"$\chi$")
    radwaves.append(radwav)
rs=np.linspace(100,10000,100)
EnergiesDelta=np.linspace(.1,4,1000)
for L in Ls:
    deltasen=[]
    sindeltasen=[]
    Smats=[]
    rmats=[]
    rmatsfromsmats=[]
    crosssections=[]
    i=0
    for En in EnergiesDelta:
        i=i+1
        radwav=chi(r,En,L,init,ain)
        radwav.SetDiffEq()
        radwav.SolveDifeq()
        radwav.SetRMatrix()
        radwav.SetSMatrix()
        radwav.SetDelta()
        radwav.SetSinDelta()
        radwav.SetCrossSection()
        deltasen.append(radwav.delta)
        sindeltasen.append(radwav.sindelta)
        Smats.append(radwav.SMat)
        rmats.append(np.abs(radwav.Rmat))
        crosssections.append(radwav.crosssection)
        rmatsfromsmats.append(1/radwav.rvalues[ain]*(HenkelMinus(radwav.GetK()*radwav.rvalues[ain],radwav.L)-radwav.SMat*HenkelPlus(radwav.GetK()*radwav.rvalues[-1],radwav.L))/((HenkelMinusPrime(radwav.GetK()*radwav.rvalues[-1],radwav.L,radwav.GetK())-radwav.SMat*HenkelPlusPrime(radwav.GetK()*radwav.rvalues[-1],radwav.L,radwav.GetK()))))
    twodplot(EnergiesDelta,makeanglecontinuous(deltasen),r" $\delta$ vs Energy for  L "+str(L),"E(MeV)", r"$\delta$")
    twodplot(EnergiesDelta,crosssections," Cross section vs Energy for  L "+str(L),"E (MeV)", " Cross Section (barns)")
    if(L==2):
        print(findpeakinfo(EnergiesDelta[250:],crosssections))
plt.show()

