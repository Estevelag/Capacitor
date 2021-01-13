made by: Andres Yesid Moreno and Esteban Velásquez

from numpy import asarray
from numpy import savetxt
import matplotlib.pyplot as plt
import numpy as np
import os.path
from PIL import Image
from copy import copy, deepcopy


#--------------------------------------------------------USER SETUP---------------------------------------------------
BoundaryType=1 #Choose: 0 for Neumann boundary conditions or 1 for Dirichlet boundary conditions
GradType=1 #Choose: 0 for manual gradient calculation or 1 for automathic gradient calculation
SF=5 #Scaling Factor (matrix elements per mm)
d=1 #separation between plates (mm)
D=260 #Diameter of Plates (mm)
UPV= 50 #upper plate voltage (V)
LPV= -50 #Lower plate Voltage (V)
Ethr=0.001 #Error threshold admitted for error between iterations


#------------------------------------------------------BACKGROUND SETUP------------------------------------------------
UDSetup= (round(SF*(d+2*d/0.5)),round(SF*(D+2*D/2))) #Universe Dimensions Set up (SF*mm)
UDSetupP=((d+2*d/0.5),(D+2*D/2)) #Copy of original UDSetup but without rounding
Universe= np.zeros(UDSetup, dtype=float) #Universe Matrix
UniverseTest=np.ones(UDSetup, dtype=float) #UniverseTest Matrix (used in order to find plate limits)
Bordex=len(Universe[0]) #Upper Universe-bound in x
Bordey=len(Universe) #Upper Universe-bound in y
Spx=round(Bordex/2) #Start position in x
Spy=round(Bordey/2) #Start position in y
iteration=0 #Count of iterations undergone by the simulation
c=0 #Logical variable
yymin=0
yymax=0
xxmin1=0
xxmin2=0
xxmax1=0
xxmax2=0


#------------------------------------------------------------UNIVERSE INITIALIZATION----------------------------------------

for i in range(SF*round(D/2)):
    Universe[Spy+round(SF*d/2),Spx+i]=UPV   #Here the electric potential of both plates is initialized
    Universe[Spy+round(SF*d/2),Spx-i]=UPV
    Universe[Spy-round(SF*d/2),Spx+i]=LPV
    Universe[Spy-round(SF*d/2),Spx-i]=LPV    
    
    
    
#-----------------------------------------------------------PLATE LIMITS CALCULATION--------- EG.(Y_UPPER,Y_LOWER, XMAX, XMIN)---------------------------------------------------    



for i in range(SF*round(D/2)):
    UniverseTest[Spy+round(SF*d/2),Spx+i]=0
    UniverseTest[Spy+round(SF*d/2),Spx-i]=0
    UniverseTest[Spy-round(SF*d/2),Spx+i]=0
    UniverseTest[Spy-round(SF*d/2),Spx-i]=0
    
cc=0
for y in range(Bordey):
    for x in range(Bordex):
        
        if (UniverseTest[y][x]==0 and cc==0):               #This section of code goes over the Universe matrix, trying to figure out the... 
            cc=1                                           #...exact discrete limits of the plates (constructed above) inside said matrix. 
            yymin=y
            xxmin1=x
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









#---------------------------------------------------------------ELECTRIC POTENTIAL FIELD CALCULATION (UNIVERSE ITERATION)-------------------------------------------------------------------




if (BoundaryType==0): #NEUMANN CONDITIONS --------------------------------------------

    while c!=1:
        
        for y in range(Bordey):
            for x in range(Bordex):
                
                if ( y>0 and y<(Bordey-1) and x>0 and x<(Bordex-1) ): 
                    
                    if  (y==Y_lower and x>=Xmin and x<=Xmax ) or (y==Y_upper and x>=Xmin and x<=Xmax ):
                     
                     UniverseC[y][x]=Universe[y][x]
                     
                    else: UniverseC[y][x]=0.25*(Universe[y][x-1]+Universe[y][x+1]+Universe[y-1][x]+Universe[y+1][x])
                     
                     
                elif (y==0 and x==0):
                     
                     UniverseC[y][x]=0.5*(Universe[y+1][x]+Universe[y][x+1])
                     
                     
                elif (y==(Bordey-1) and x==0):
                     
                     UniverseC[y][x]=0.5*(Universe[y-1][x]+Universe[y][x+1])           #For the Neumann boundary conditions, the value of the electric potential of a cell at the boundary is given by the average valueof its neighbours
                                                                                       #Aside from that, the algorithym for computing the rest of the electric potential is the same.   
                                                                                       #Two matrices of electric potential (Universe and UniverseC) are used together in order to aid the computing process.
                elif (y==0 and x==(Bordex-1)):
                     
                     UniverseC[y][x]=0.5*(Universe[y+1][x]+Universe[y][x-1])
                     
                     
                elif (y==(Bordey-1) and x==(Bordex-1)):
                     
                     UniverseC[y][x]=0.5*(Universe[y-1][x]+Universe[y][x-1])
                     
                     
                elif (y==0 and x>0 and x<(Bordex-1)):
                     
                     UniverseC[y][x]=0.3333*(Universe[y][x+1]+Universe[y][x-1]+Universe[y+1][x])           #Here the matrix for UniverseC is updated using the values in the matrix Universe
                     
                     
                elif (y==(Bordey-1) and x>0 and x<(Bordex-1)):
                     
                     UniverseC[y][x]=0.3333*(Universe[y][x+1]+Universe[y][x-1]+Universe[y-1][x])
                     
                     
                elif (x==0 and y>0 and y<(Bordey-1)):
                     
                     UniverseC[y][x]=0.3333*(Universe[y+1][x]+Universe[y-1][x]+Universe[y][x+1])  
                     
                     
                elif (x==(Bordex-1) and y>0 and y<(Bordey-1)):
                     
                     UniverseC[y][x]=0.3333*(Universe[y+1][x]+Universe[y-1][x]+Universe[y][x-1])  
                     
        iteration=iteration+1            
        
        for y in range(Bordey):
            for x in range(Bordex):
                
                if ( y>0 and y<(Bordey-1) and x>0 and x<(Bordex-1) ): 
                    
                    if  (y==Y_lower and x>=Xmin and x<=Xmax ) or (y==Y_upper and x>=Xmin and x<=Xmax ):
                     
                     Universe[y][x]=UniverseC[y][x]
                     
                    else: Universe[y][x]=0.25*(UniverseC[y][x-1]+UniverseC[y][x+1]+UniverseC[y-1][x]+UniverseC[y+1][x])
                     
                     
                elif (y==0 and x==0):
                     
                     Universe[y][x]=0.5*(UniverseC[y+1][x]+UniverseC[y][x+1])
                     
                     
                elif (y==(Bordey-1) and x==0):
                     
                     Universe[y][x]=0.5*(UniverseC[y-1][x]+UniverseC[y][x+1])
                     
                     
                elif (y==0 and x==(Bordex-1)):
                     
                     Universe[y][x]=0.5*(UniverseC[y+1][x]+UniverseC[y][x-1])                                   #Here the matrix for Universe is updated using the values in the matrix UniverseC
                     
                     
                elif (y==(Bordey-1) and x==(Bordex-1)):
                     
                     Universe[y][x]=0.5*(UniverseC[y-1][x]+UniverseC[y][x-1])
                     
                     
                elif (y==0 and x>0 and x<(Bordex-1)):
                     
                     Universe[y][x]=0.3333*(UniverseC[y][x+1]+UniverseC[y][x-1]+UniverseC[y+1][x])
                     
                     
                elif (y==(Bordey-1) and x>0 and x<(Bordex-1)):
                     
                     Universe[y][x]=0.3333*(UniverseC[y][x+1]+UniverseC[y][x-1]+UniverseC[y-1][x])
                     
                     
                elif (x==0 and y>0 and y<(Bordey-1)):
                     
                     Universe[y][x]=0.3333*(UniverseC[y+1][x]+UniverseC[y-1][x]+UniverseC[y][x+1])  
                     
                     
                elif (x==(Bordex-1) and y>0 and y<(Bordey-1)):
                     
                     Universe[y][x]=0.3333*(UniverseC[y+1][x]+UniverseC[y-1][x]+UniverseC[y][x-1])   
                     
        iteration=iteration+1 
        
        if (iteration>Bordey and ((iteration%50)==0)):
            errormax1=0
            Error=abs(Universe-UniverseC)
            for i in range(len(Error)):
               errormax2= max(Error[:][i])              #This section of code computes the convergence error between iterations. When the error threshold is reached, the iteration process stops
               errormax1=max(errormax1,errormax2) 
            print(errormax1)
               
            if errormax1<Ethr:
             c=1;
              
        print(iteration)           


elif (BoundaryType==1): #DIRICHLET CONDITIONS ----------------------------------------------------
    
        while c!=1:
                                                                                     #For the Dirichlet boundary conditions, the value of the electric potential of a cell at the boundary is constant through time and equal to zero      
            for y in range(Bordey):                                                  #Aside from that, the algorithym for computing the rest of the electric potential is the same.
                for x in range(Bordex):
                    
                    if ( y>0 and y<(Bordey-1) and x>0 and x<(Bordex-1) ): 
                        
                        if  (y==Y_lower and x>=Xmin and x<=Xmax ) or (y==Y_upper and x>=Xmin and x<=Xmax ):
                         
                         UniverseC[y][x]=Universe[y][x]
                         
                        else: UniverseC[y][x]=0.25*(Universe[y][x-1]+Universe[y][x+1]+Universe[y-1][x]+Universe[y+1][x])
                         
            iteration=iteration+1            
            
            for y in range(Bordey):
                for x in range(Bordex):
                    
                    if ( y>0 and y<(Bordey-1) and x>0 and x<(Bordex-1) ): 
                        
                        if  (y==Y_lower and x>=Xmin and x<=Xmax ) or (y==Y_upper and x>=Xmin and x<=Xmax ):
                         
                         Universe[y][x]=UniverseC[y][x]
                         
                        else: Universe[y][x]=0.25*(UniverseC[y][x-1]+UniverseC[y][x+1]+UniverseC[y-1][x]+UniverseC[y+1][x])
                         
            iteration=iteration+1 
            
            if (iteration>Bordey and ((iteration%50)==0)):
                errormax1=0
                Error=abs(Universe-UniverseC)
                for i in range(len(Error)):
                   errormax2= max(Error[:][i])
                   errormax1=max(errormax1,errormax2) 
                print(errormax1)
                   
                if errormax1<Ethr:
                 c=1;
                  
            print(iteration)         


#---------------------------------------------------UNIVERSE MATRIX SAVE-----------------------------------------------------------

StringSave='Universe_2mm_100V_0.05%_Dirichlet_V1.csv'

if ( (d!=Ethr*100) and (d!=abs(UPV-LPV)) and (Ethr*100!=abs(UPV-LPV)) ):

    StringSave=StringSave.replace('2',str(d))
    StringSave=StringSave.replace('0.05',str(Ethr*100))
    StringSave=StringSave.replace('100',str(abs(UPV-LPV)))          #This section of code just automates the saving of Electric Potential matrices...
                                                                   #...according to some previously-given nomenclature criteria 
    if (BoundaryType==1):
        StringSave=StringSave.replace('Dirichlet','Dirichlet')
    elif (BoundaryType==0):
        StringSave=StringSave.replace('Dirichlet','Neumann')
        
    String1=StringSave
    ex=1
    String2='V1'
    while os.path.exists(String1)==True:
        ex=ex+1;
        String3=String2.replace('1',str(ex))
        String3p=String2.replace('1',str(ex-1))
        String1=String1.replace(String3p,String3)
        
    StringSave=String1    
else:
    print('An unforseen situation has occurred!')
    StringSave=input('Please type manually the name with which you want to save the Electric Potential (Universe) matrix:')
        
    
savetxt(StringSave, Universe, delimiter=',')


StingSetUp='UDSetup--M'
StringSetUpP='UDSetupP--M'
StringSetUp=StingSetUp.replace('M',StringSave)
StringSetUpP=StringSetUpP.replace('M',StringSave)
savetxt(StringSetUp,UDSetup, delimiter=',')           #Saves data of the calculated electric field for postprocessing
savetxt(StringSetUpP,UDSetupP, delimiter=',')




