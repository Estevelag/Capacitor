from numpy import asarray
from numpy import savetxt
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import numpy as np
from PIL import Image
from copy import copy, deepcopy
from numpy import loadtxt
from scipy import interpolate
import pylab as py
import math
import scipy.constants as const
#--------------------------------------------------------USER SETUP--------------------------------------------------------------
#StringLoad=input('Please enter the file name of the Electric Potential data that you would like to use:')
StringLoad='Universe_2mm_100V_0.1%_Dirichlet_V1.csv'
StingSetUp='UDSetup--M'
StringSetUpP='UDSetupP--M'
StringSetUp=StingSetUp.replace('M',StringLoad)
StringSetUpP=StringSetUpP.replace('M',StringLoad)

BoundaryType=1 #Choose: 0 for Neumann boundary conditions or 1 for Dirichet boundary conditions
GradType=1 #Choose: 0 for manual gradient calculation or 1 for automathic gradient calculation
SF=5 #Scaling Factor
d=2 #separation between plates [mm]
D=260 #Diameter of Plates [mm]
UPV=50 #upper plate voltage [V]
LPV=-50 #Lower plate Voltage [V]

#-----------------------------------------------------BACKGROUND SETUP------------------------------------------------------------------
UDSetup= loadtxt(StringSetUp, delimiter=',')
UDSetupP=loadtxt(StringSetUpP, delimiter=',')
UDSetup=(int(UDSetup[0]),int(UDSetup[1]))
UDSetupP=(UDSetupP[0],UDSetupP[1])
Universe= np.zeros(UDSetup, dtype=float) #Universe Matrix
UniverseTest=np.ones(UDSetup, dtype=float) #UniverseTest Matrix (used in order to find plate limits)
Bordex=len(Universe[0]) #Upper Universe-bound in x
Bordey=len(Universe) #Upper Universe-bound in y
Spx=round(Bordex/2) #Start position in x
Spy=round(Bordey/2) #Start position in y
iteration=0 #Iteration counter
c=0 #Logical Variable
Eo=8.8541878128*10**(-12) #Vacuum permittivity [F/m]
yymin=0
yymax=0
xxmin1=0
xxmin2=0
xxmax1=0
xxmax2=0

    
#-----------------------------------------------------------PLATE LIMITS SET UP------------------------------------------------------------    

for i in range(SF*round(D/2)):
    UniverseTest[Spy+round(SF*d/2),Spx+i]=0
    UniverseTest[Spy+round(SF*d/2),Spx-i]=0
    UniverseTest[Spy-round(SF*d/2),Spx+i]=0
    UniverseTest[Spy-round(SF*d/2),Spx-i]=0
    
cc=0
for y in range(Bordey):
    for x in range(Bordex):
        
        if (UniverseTest[y][x]==0 and cc==0):
            cc=1
            yymin=y                                                              #This section of code goes over the Universe matrix, trying to figure out the...
            xxmin1=x                                                             #...exact discrete limits of the plates (constructed above) inside said matrix. 
        elif (y==yymin and x>=xxmin1 and UniverseTest[y][x]==1 and cc==1):
            xxmax1=x-1
            cc=2
        elif (UniverseTest[y][x]==0 and y!=yymin and cc==2):
            yymax=y
            xxmin2=x
            cc=3
        elif (y==yymax and x>=xxmin2 and UniverseTest[y][x]==1 and cc==3):
            xxmax2=x-1
            cc=4

if(xxmin1==xxmin2 and xxmax1==xxmax2):
    Y_lower= yymin
    Y_upper= yymax
    Xmin= xxmin1
    Xmax= xxmax1
    
UniverseC=deepcopy(Universe)


#-----------------------------------------------------------------------UNIVERSE LOAD------------------------------------------------

Universe = loadtxt(StringLoad, delimiter=',')
#Loads precalculated electric potential field data

#---------------------------------------------------------------------CALCULATION OF ELECTRIC FIELD- AND CAPACITANCE--------------------------


if (GradType==0): #MANUAL GRADIENT AND CAPACITANCE--------------

    UniverseEx=np.zeros(UDSetup, dtype=float) 
    UniverseEy=np.zeros(UDSetup, dtype=float)
                                                      #Se realiza la diferencia centrada
    for y in range(Bordey):
        for x in range(Bordex):
            
            if x>0 and x<(Bordex-1) :
                UniverseEx[y][x]=(Universe[y][x+1]-Universe[y][x-1])/(2*(1/SF)*10**(-3)) 
            elif x==0:
                 UniverseEx[y][x]=(Universe[y][x+1]-0)/(2*(1/SF)*10**(-3))   
            elif x==(Bordex-1):
                 UniverseEx[y][x]=(0-Universe[y][x-1])/(2*(1/SF)*10**(-3))  
                 
            
    for y in range(Bordey):
        for x in range(Bordex):
            
            if y>0 and y<(Bordey-1):
                UniverseEy[y][x]=(Universe[y+1][x]-Universe[y-1][x])/(2*(1/SF)*10**(-3)) 
            elif y==0:
                 UniverseEy[y][x]=(Universe[y+1][x]-0)/(2*(1/SF)*10**(-3)) 
            elif y==(Bordey-1):
                 UniverseEy[y][x]=(0-Universe[y-1][x])/(2*(1/SF)*10**(-3))     
                
    UniverseEx=-UniverseEx
    UniverseEy=-UniverseEy     
                           
    z=0
    Charge1=0
    while (z!=Xmax):
        Charge1= Charge1 + abs(UniverseEy[Y_upper-1][Xmin+z])*(1/SF)*(10**(-3))*Eo
        z=z+1
    Charge1=Charge1+(abs(UniverseEx[Y_upper][Xmin-1]) + abs(UniverseEx[Y_upper][Xmax+1]))*(1/SF)*(10**(-3))*Eo
    Charge1=(Charge1/2)*2*np.pi*((D/2)*10**(-3))**2
    
    Capacitance1= (Charge1/(abs(UPV-LPV)))*10**(9) #Capacitance in nanoFarads
    print('Capacitance1 [nF]=') 
    print(Capacitance1)
    
    #this is another way to calculate the capacitance using geometrydifferential areas
    Qt=0
    Areatot=0
    for i in range(Xmin,round((Xmax-Xmin)/2)+Xmin+1):#first we take the first half of the plate and by simmetry this is half the area and the electric field behaves equally
        R=(round((Xmax-Xmin)/2))+1
        D=((round((Xmax-Xmin)/2))+1-(i-Xmin+1))# this is the distance from the center of the plate that is the loop currently in spaces
        Theta=(180/math.pi)*math.acos(D/R)#This is the first angle in which we are in
        Theta2=(180/math.pi)*math.acos((D+1)/R)# This is the next angle in the next value
        Area=4*(((R**2*math.pi*(Theta*2)/360)-(R*D*math.sin(Theta*math.pi/180)))-(((R)**2*math.pi*(Theta2*2)/360)-((D+1)*R*math.sin(Theta2*math.pi/180))))
        #The area denotes the area in 0.1mm units of a slice of circle with the width of 0.2mm
        Qt=Qt+abs((UniverseEy[Y_lower+1,i]**2+UniverseEx[Y_lower+1,i]**2)**0.5)*Area/(10000**2)#This aproximates that the field in the perpendicular axe to the condensator is equal and the 10000 is to make it in SI
        Areatot=Areatot+Area# This defines the total area utilized
    #Now we calculate the capacitance with this aproximation taking into account the border
    Capacitance3=10*Qt*const.epsilon_0/100
    Areasobred=math.pi*0.13**2/(0.001*d)
    Err=100*(Capacitance3-Areasobred*const.epsilon_0)/(Areasobred*const.epsilon_0)
    Areatot=2*Areatot/(10000**2)#This is to verify if it really made the whole circle and the factor is to convert to 0.1mm in area ab because is half the area
    print('Capacitance3 [nF]=') 
    print(Capacitance3*10**(9))
    
    
   
elif (GradType==1): #AUTOMATIC GRADIENT AND CAPACITANCE

    
    Vgrad=np.gradient(Universe)
    UniverseEx=-Vgrad[1]*(10**3)
    UniverseEy=-Vgrad[0]*(10**3)
    
#    z=0
#    Capacitance2=0
#    while (z!=Xmax):
#        Capacitance2= Capacitance2 + abs(UniverseEy[Y_upper-1][Xmin+z])*(1/SF)*(10**(-3))*Eo
#        z=z+1
#        Capacitance2=Capacitance2 + (abs(UniverseEx[Y_upper][Xmin-1]) + abs(UniverseEx[Y_upper][Xmax+1]))*(1/SF)*(10**(-3))*Eo
#        Capacitance2=(Capacitance2/2)*2*np.pi*(((D/2)*10**(-3))**2)
#    Capacitance2= (Capacitance2/(abs(UPV-LPV)))*10**9 #Capacitance in nanoFarads
#    print('Capacitance2 [nF]=') 
#    print(Capacitance2)
    
    z=0
    Capacitance2=0
    while (z!=Xmax):
        Capacitance2= Capacitance2 + (((UniverseEx[Y_upper-1][Xmin+z])**2 + (UniverseEy[Y_upper-1][Xmin+z])**2)**0.5)*(1/SF)*(10**(-3))*Eo
        z=z+1
    
    Capacitance2= (Capacitance2/(abs(UPV-LPV)))*10**(9) #Capacitance in nanoFarads
    print('Capacitance2 [nF]=') 
    print(Capacitance2)
    
    #this is another way to calculate the capacitance using geometrical differential areas
    Qt=0
    Areatot=0
    for i in range(Xmin,round((Xmax-Xmin)/2)+Xmin+1):#first we take the first half of the plate and by simmetry this is half the area and the electric field behaves equally
        R=(round((Xmax-Xmin)/2))+1
        D=((round((Xmax-Xmin)/2))+1-(i-Xmin+1))# this is the distance from the center of the plate that is the loop currently in spaces
        Theta=(180/math.pi)*math.acos(D/R)#This is the first angle in which we are in
        Theta2=(180/math.pi)*math.acos((D+1)/R)# This is the next angle in the next value
        Area=4*(((R**2*math.pi*(Theta*2)/360)-(R*D*math.sin(Theta*math.pi/180)))-(((R)**2*math.pi*(Theta2*2)/360)-((D+1)*R*math.sin(Theta2*math.pi/180))))
        #The area denotes the area in 0.1mm units of a slice of circle with the width of 0.2mm
        Qt=Qt+abs((UniverseEy[Y_lower+1,i]**2+UniverseEx[Y_lower+1,i]**2)**0.5)*Area/(10000**2)#This aproximates that the field in the perpendicular axe to the condensator is equal and the 10000 is to make it in SI
        Areatot=Areatot+Area# This defines the total area utilized
    #Now we calculate the capacitance with this aproximation taking into account the border
    Capacitance3=10*Qt*const.epsilon_0/100
    Areasobred=math.pi*0.13**2/(0.001*d)
    Err=100*(Capacitance3-Areasobred*const.epsilon_0)/(Areasobred*const.epsilon_0)
    Areatot=2*Areatot/(10000**2)#This is to verify if it really made the whole circle and the factor is to convert to 0.1mm in area ab because is half the area
    print('Capacitance3 [nF]=') 
    print(Capacitance3*10**(9))
    



  
#------------------------------------------------------DATA VISUALIZATION----------------------------------------------
Emagnitude=np.power(UniverseEx,2)+np.power(UniverseEy,2)
Emagnitude=np.power(Emagnitude,0.5)
      
UDSetupP=(int(UDSetupP[0]),int(UDSetupP[1]))
# Grid of x, y points
nx, ny = SF*UDSetupP[1], SF*UDSetupP[0]
x = np.linspace(0, UDSetupP[1], nx)
y = np.linspace(0, UDSetupP[0], ny)
X, Y = np.meshgrid(x, y)


fx = interpolate.interp2d(x, y, UniverseEx, kind='cubic')
fy = interpolate.interp2d(x, y, UniverseEy, kind='cubic')

#EQUIPOTENTIAL LINES VISUALIZATION    
breaks=np.linspace(LPV, UPV ,21)
plt.figure(figsize=(9,8))
CS1 = plt.contourf(X, Y, Universe, 
breaks,
cmap='jet' )
plt.colorbar(ticks=breaks, orientation= 'vertical')
plt.xlabel('Horizontal Distance [mm]')    
plt.ylabel('Vertical Distance [mm]')
plt.title('Equipotential lines of a 2D parallel plate capacitor [V].')         


#ELECTRIC FIELD LINES VISUALIZATION
fig1=plt.figure(figsize=(9,8))      
ax = fig1.add_subplot(1, 1, 1) 
rect1 = plt.Rectangle((Xmin/SF, Y_upper/SF), D, 1/SF, color='k', alpha=0.5)   
rect2 = plt.Rectangle((Xmin/SF, Y_lower/SF), D, 1/SF, color='k', alpha=0.5)     
plt.streamplot(X, Y, fx(x,y), fy(x,y), color=Emagnitude, density=[1,10], cmap=plt.cm.jet) 
plt.xlabel('Horizontal Distance [mm]')    
plt.ylabel('Vertical Distance [mm]')
plt.title('Electric field lines of  a 2D parallel plate capacitor [V/m]')     
plt.colorbar(orientation= 'vertical')
ax.add_patch(rect1)
ax.add_patch(rect2)

#ELECTRIC FIELD INTENSITY VISUALIZATION
breaks2=np.linspace(0, np.amax(Emagnitude), 25)
plt.figure(figsize=(9,8))
CS1 = plt.contourf(X, Y, Emagnitude, breaks2,
cmap='jet' )
plt.colorbar(ticks=breaks2, orientation= 'vertical')
plt.xlabel('Horizontal Distance [mm]')    
plt.ylabel('Vertical Distance [mm]')
plt.title('Electric field intensity of  a 2D parallel plate capacitor [V/m]') 





   
 



























#    z=0
#    Capacitance2=0
#    while (z!=Xmax):
#        Capacitance2= Capacitance2 + (((UniverseEx[Y_upper-1][Xmin+z])**2 + (UniverseEy[Y_upper-1][Xmin+z])**2)**0.5)*(1/SF)*(10**(-3))*Eo
#        z=z+1
#    
#    Capacitance2= (Capacitance2/(abs(UPV-LPV)))*10**(9) #Capacitance in nanoFarads
#    print('Capacitance2 [nF]=') 
#    print(Capacitance2)
