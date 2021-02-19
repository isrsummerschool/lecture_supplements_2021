"""

Routines to provide collision frequency coefficients, from 
Schunk and Nagy [2000]

M Nicolls - 2013

"""

import scipy

# getBst
def getBst(s,t):
    '''
    
    Ion-coulomb collision frequency coefficients from Schunk and Nagy [2000] Table 4.3.
    Collision frequencies can then be calculated using Eq. 4.143 in SN00.
    
    Inputs:
        s,t - in amu
        (1=H+,4=He+,12=C+,14=N+,16=O+,29=CO+,28=N2+,30=NO+,32=O2+,44=CO2+)
    
    Note for CO+, need to set s or t to 29 amu (because of duplicate N2+ value) 
    '''
        
    Amu2Ion = {'1':'H+','4':'He+','12':'C+','14':'N+','16':'O+','29':'CO+','28':'N2+','30':'NO+','32':'O2+','44':'CO2+'}

    Bst = {
        'H+':   {'H+':0.90, 'He+':1.14, 'C+':1.22,'N+':1.23,'O+':1.23,'CO+':1.25,'N2+':1.25,'NO+':1.25,'O2+':1.25,'CO2+':1.26},
        'He+':  {'H+':0.28, 'He+':0.45, 'C+':0.55,'N+':0.56,'O+':0.57,'CO+':0.59,'N2+':0.59,'NO+':0.60,'O2+':0.60,'CO2+':0.61},
        'C+':   {'H+':0.102,'He+':0.18, 'C+':0.26,'N+':0.27,'O+':0.28,'CO+':0.31,'N2+':0.31,'NO+':0.31,'O2+':0.31,'CO2+':0.32},
        'N+':   {'H+':0.088,'He+':0.16, 'C+':0.23,'N+':0.24,'O+':0.25,'CO+':0.28,'N2+':0.28,'NO+':0.28,'O2+':0.28,'CO2+':0.30},
        'O+':   {'H+':0.077,'He+':0.14, 'C+':0.21,'N+':0.22,'O+':0.22,'CO+':0.25,'N2+':0.25,'NO+':0.26,'O2+':0.26,'CO2+':0.27},
        'CO+':  {'H+':0.045,'He+':0.085,'C+':0.13,'N+':0.14,'O+':0.15,'CO+':0.17,'N2+':0.17,'NO+':0.17,'O2+':0.18,'CO2+':0.19},
        'N2+':  {'H+':0.045,'He+':0.085,'C+':0.13,'N+':0.14,'O+':0.15,'CO+':0.17,'N2+':0.17,'NO+':0.17,'O2+':0.18,'CO2+':0.19},
        'NO+':  {'H+':0.042,'He+':0.080,'C+':0.12,'N+':0.13,'O+':0.14,'CO+':0.16,'N2+':0.16,'NO+':0.16,'O2+':0.17,'CO2+':0.18},
        'O2+':  {'H+':0.039,'He+':0.075,'C+':0.12,'N+':0.12,'O+':0.13,'CO+':0.15,'N2+':0.15,'NO+':0.16,'O2+':0.16,'CO2+':0.17},
        'CO2+': {'H+':0.029,'He+':0.055,'C+':0.09,'N+':0.09,'O+':0.10,'CO+':0.12,'N2+':0.12,'NO+':0.12,'O2+':0.12,'CO2+':0.14},
        }
        
    try:
        s=int(s)
        t=int(t)
        n1=Amu2Ion[str(s)]
        n2=Amu2Ion[str(t)]
        return Bst[n1][n2]
    except:
        return 0.0

# getCin
def getCin(i,n,Ti=1000.0,Tn=1000.0):
    '''
    Non-resonant and resonant collision frequency coefficients (Cin x 10^10) from SN00 Table 4.4 and 4.5
    Collision frequencies can then be calculated as vin = Cin x 1e-10 x Nn where Nn is the neutral density in cm^-3
    
    Inputs:
        i - ion mass in amu (1=H+,4=He+,12=C+,14=N+,16=O+,29=CO+,28=N2+,30=NO+,32=O2+,44=CO2+)
        n - neutral mass in amu (1=H,4=He,14=N,16=O,29=CO,28=N2,32=O2,44=CO2)
        Ti - ion temp (for resonant collisions)
        Tn - neutral temp (for resonant collisions)

    Note for CO+, need to set i or n to 29 amu (because of duplicate N2+ value) 
    
    '''
    
    Amu2Ion = {'1':'H+','4':'He+','12':'C+','14':'N+','16':'O+','29':'CO+','28':'N2+','30':'NO+','32':'O2+','44':'CO2+'}
    Amu2Ntrl = {'1':'H','4':'He','14':'N','16':'O','29':'CO','28':'N2','32':'O2','44':'CO2'}

    # resonant
    Tr = (Ti+Tn)/2.0
    HpH = 2.65*Tr**0.5*(1.0-0.083*scipy.log10(Tr))**2.0
    HpO = 0.661*Ti**0.5*(1.0-0.047*scipy.log10(Ti))**2.0
    HepHe = 0.873*Tr**0.5*(1.0-0.093*scipy.log10(Tr))**2.0
    NpN = 0.383*Tr**0.5*(1.0-0.063*scipy.log10(Tr))**2.0
    OpH = 0.661*Ti**0.5*(1.0-0.047*scipy.log10(Ti))**2.0
    OpO = 0.367*Tr**0.5*(1.0-0.064*scipy.log10(Tr))**2.0
    COpCO = 0.342*Tr**0.5*(1.0-0.085*scipy.log10(Tr))**2.0
    N2pN2 = 0.514*Tr**0.5*(1.0-0.073*scipy.log10(Tr))**2.0
    O2pO2 = 0.259*Tr**0.5*(1.0-0.063*scipy.log10(Tr))**2.0
    CO2pCO2 = 0.285*Tr**0.5*(1.0-0.083*scipy.log10(Tr))**2.0
    
    Cin = {
        'H+':   {'H':HpH, 'He':10.6, 'N':26.1,'O':HpO, 'CO':35.6, 'N2':33.6, 'O2':32.0, 'CO2':41.4},
        'He+':  {'H':4.71,'He':HepHe,'N':11.9,'O':10.1,'CO':16.9, 'N2':16.0, 'O2':15.3, 'CO2':20.0},
        'C+':   {'H':1.69,'He':1.71, 'N':5.73,'O':4.94,'CO':8.74, 'N2':8.26, 'O2':8.01, 'CO2':10.7},
        'N+':   {'H':1.45,'He':1.49, 'N':NpN, 'O':4.42,'CO':7.90, 'N2':7.47, 'O2':7.25, 'CO2':9.73},
        'O+':   {'H':OpH, 'He':1.32, 'N':4.62,'O':OpO, 'CO':7.22, 'N2':6.82, 'O2':6.64, 'CO2':8.95},
        'CO+':  {'H':0.74,'He':0.79, 'N':2.95,'O':2.58,'CO':COpCO,'N2':4.24, 'O2':4.49, 'CO2':6.18},
        'N2+':  {'H':0.74,'He':0.79, 'N':2.95,'O':2.58,'CO':4.84, 'N2':N2pN2,'O2':4.49, 'CO2':6.18},
        'NO+':  {'H':0.69,'He':0.74, 'N':2.79,'O':2.44,'CO':4.59, 'N2':4.34, 'O2':4.27, 'CO2':5.89},
        'O2+':  {'H':0.65,'He':0.70, 'N':2.64,'O':2.31,'CO':4.37, 'N2':4.13, 'O2':O2pO2,'CO2':5.63},
        'CO2+': {'H':0.47,'He':0.51, 'N':2.00,'O':1.76,'CO':3.40, 'N2':3.22, 'O2':3.18, 'CO2':CO2pCO2},
        }

    try:
        i=int(i)
        n=int(n)
        ion=Amu2Ion[str(i)]
        ntrl=Amu2Ntrl[str(n)]        
        return Cin[ion][ntrl]
    except:
        return 0.0
        
# getCen
def getCen(n,Te=1000.0):

    '''
    Electron-neutral collision frequency coefficients (Cen x 10^10) from SN00 Table 4.6
    Collision frequencies can then be calculated as ven = Cen x 1e-10 x Nn where Nn is the neutral density in cm^-3
    
    Inputs:
        n - neutral mass in amu (1=H,4=He,16=O,29=CO,28=N2,32=O2,44=CO2)
        Te - electron temp 

    Note for CO+, need to set e to 29 amu (because of duplicate N2+ value) 
    
    '''

    Amu2Ntrl = {'1':'H','4':'He','16':'O','29':'CO','28':'N2','32':'O2','44':'CO2'}

    Cen = {
        'N2':0.233*(1.0-1.21e-4*Te)*Te, 
        'O2':1.82*(1.0+3.6e-2*Te**0.5)*Te**0.5, 
        'O':0.89*(1.0+5.7e-4*Te)*Te**0.5, 
        'He':0.46*Te**0.5,
        'H':45.0*(1.0-1.35e-4*Te)*Te**0.5, 
        'CO':0.234*(Te+165.0), 
        'CO2':368.0*(1.0+4.1e-11*scipy.absolute(4500.0-Te)**2.93),
    }

    try:
        n=int(n)
        ntrl=Amu2Ntrl[str(n)]        
        return Cen[ntrl]
    except:
        return 0.0

    return