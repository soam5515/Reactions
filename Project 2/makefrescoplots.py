import numpy as np
import matplotlib.pyplot as plt
import os


def GetCrossSectionData(f):
    arr=np.genfromtxt(f,skiprows=14,skip_footer=1)
    return arr

def twodplot(x,y,title,xaxis,yaxis):
    fig=plt.figure()
    ax2=fig.add_subplot(111)
    ax2.plot(x,y)
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
#l="Ni58coulombpointlikepoint1MeV"
#h="pointlikeNi58CoulombEn50MeV"
#hc="Ni58coulombEn50MeV"
#lc="Ni58coulombpoint1MeV"

#files=[]
di='C:\\Users\\Charles\\Documents\\frescoreactions\\reactions\\'

Neu50=GetCrossSectionData(di+"Ni58Neutron50MeV\\fort.16")
Neu5=GetCrossSectionData(di+"Ni58Neutron5MeV\\fort.16")

twodplot(Neu50[:,0],Neu50[:,1],"Neutron-Nickel 58 Cross Section for E=50 MeV","Angle (Degrees)", "Differential Cross Section(relative to Rutherford")
twodplot(Neu5[:,0],Neu5[:,1],"Neutron-Nickel 58 Cross Section for E=5 MeV","Angle (Degrees)", "Differential Cross Section(relative to Rutherford")



#for folder in os.listdir(di):
#    files.append(di+folder+"\\fort.16")
#print(os.listdir(di))
#f=files[1]
#a=np.loadtxt(files[1],comments='@',usecols=(0,1))
#a=np.genfromtxt(f,comments='@',skiprows=10,skip_footer=1)
# lowenpoint=GetCrossSectionData(di+l+"\\fort.16")
# twodplot(lowenpoint[:,0],lowenpoint[:,1],"Pointlike Coulomb Scattering p-Ni58, E=.1 MeV","Angle (Degrees)","Differential Cross Section (Relative to Rutherford)")
# highenpoint=GetCrossSectionData(di+h+"\\fort.16")
# twodplot(highenpoint[:,0],highenpoint[:,1],"Pointlike Coulomb Scattering p-Ni58, E=50 MeV","Angle (Degrees)","Differential Cross Section (Relative to Rutherford)")
# highc=GetCrossSectionData(di+hc+"\\fort.16")
# lowc=GetCrossSectionData(di+lc+"\\fort.16")
# twodplot(highc[:,0],highc[:,1],"Coulomb Scattering p-Ni58, E=50 MeV","Angle (Degrees)","Differential Cross Section (Relative to Rutherford)")
# twodplot(lowc[:,0],lowc[:,1],"Coulomb Scattering p-Ni58, E=.1 MeV","Angle (Degrees)","Differential Cross Section (Relative to Rutherford)")
#xs=[]
#ys=[]
#for item in a:
#    xs.append(item[0])
#    ys.append(item[1])
#plt.plot(a[:,0],a[:,1])
plt.show()
