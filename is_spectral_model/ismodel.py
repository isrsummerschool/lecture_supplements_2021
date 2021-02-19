"""

Main IS spectral class

M Nicolls - 2011-2013

"""

import sys
import scipy,scipy.constants,scipy.fftpack
import time

from sommerfeldIntegral2 import *
import collisionCoeffs
#from ..models.utils import collisionCoeffs

#import logging

class ISspec:        

    Amu2Ion = {'1':'H+','4':'He+','12':'C+','14':'N+','16':'O+','29':'CO+','28':'N2+','30':'NO+','32':'O2+','44':'CO2+'}

    def __init__(self,input,emode=[0,0,0],imode=[0,0,0],fmax=-1.0,Nfreq=100,normalize=0,override=0,czparams=(1e-6,2e5,100,10,1.0)):
        
        """
        
        Inputs:
            input: A dictionary with the following entries
                input['Nion'] - # of ions
                input['mi'] - list of ion masses in amu, e.g. [1.0,16.0,4.0], length Nion
                input['B'] - background magnetic field values in Tesla
                input['f0'] - Frequency in Hz
                input['te'] - Electron temperature in K
                input['alpha'] - aspect angle in degrees (90 for perpendicular to B)
                input['ne'] - electron density in m^-3
                input['ti'] - list of Ti values, e.g. [1000.0,1000.0,1000.0], length Nion
                input['ni'] - list of ion fractions, e.g. [0.0, 1.0, 0.0], length Nion, must add to 1
                input['ve'] - electron speed, m/s
                input['vi'] - list of ion speeds, m/s, e.g. [0.0, 0.0, 0.0], length Nion
                input['nuen'] - electron-neutral collision frequency, s^-1
                input['nuin'] - list of ion neutral collision frequencies, s^-1, e.g. [0.0, 0.0, 0.0], length Nion
            emode, imode: flags for Gordeyev integral, 0=off, 1=on
                [0] - Magnetic field
                [1] - Coulomb collisions
                [2] - Ion-neutral collisions (BGK)
            fmax: max frequency in Hz (will calculate if not provided)
            Nfreq: determines frequency resolution
            normalize: whether to normalize ACF
            override: overrides parameter consistency check
            czparams: parameters for Gordeyev integral evaluation: (tol,Nmax,maxLoops,Nstart,kmax)
                tol - tolerence criteria for convergence
                Nmax - maximum number of points in integral
                maxLoops - maximum number of loops permissible
                Nstart - starting number of points in integral
                kmax - ending value of integral bound
            
        """

        # set up logging
        #self.logger = logging.getLogger(__name__)

        # set override flag
        self.override=override
        
        # Check to make sure input parameters are valid
        # Stores variables in self.Params dictionary
        self.prepareAndCheck(input)
        
        # Copy variables into class variables
        self.emode=emode
        self.imode=imode
        self.fmax=fmax
        self.Nfreq=Nfreq
        self.normalize=normalize
        self.czparamsi=czparams 
        self.czparamse=czparams
                
        return
                
    def adjustParams(self,newParms):
        
        """
        adjustParams: Changes parameters for IS model evaluation
        
        Inputs:
            newParms: A dictionary with entries corresponding to the elements to be changed

        """

        pkeys = newParms.keys()
                
        for key in pkeys:
            try: 
                self.Params[key]
            except Exception as e:
                self.logger.error("invalid parameter %s" % key)
                self.logger.exception(e)
                raise
                          
            if key in ('ni','ti','vi','nuin','nuii'):
                self.Params[key] = scipy.array(newParms[key])
            else:
                self.Params[key] = newParms[key]
                         
        if not self.override:
            self.checkParams()
        else:
            self.logger.info('Override flag is set, so not checking')
          
        return
        
    def prepareAndCheck(self,input):
    
        """
        prepareAndCheck: 
        Creates class dictionary "self.Params" with the following keys:
                'Nion' - # of ions
                'mi' - list of ion masses in amu, e.g. [1.0,16.0,4.0], length Nion
                'B' - background magnetic field values in Tesla
                'f0' - Frequency in Hz
                'te' - Electron temperature in K
                'alpha' - aspect angle in degrees (90 for perpendicular to B)
                'ne' - electron density in m^-3
                'ti' - list of Ti values, e.g. [1000.0,1000.0,1000.0], length Nion
                'ni' - list of ion fractions, e.g. [0.0, 1.0, 0.0], length Nion, must add to 1
                've' - electron speed, m/s
                'vi' - list of ion speeds, m/s, e.g. [0.0, 0.0, 0.0], length Nion
                'nuen' - electron-neutral collision frequency, s^-1
                'nuin' - list of ion neutral collision frequencies, s^-1, e.g. [0.0, 0.0, 0.0], length Nion
                
        Will raise KeyError when input does not Nion, mi, or f0. Otherwise will set default values.
        
        Inputs:
            Initialized "input" dictionary
        """
        
        self.Params = {}
    
        try:
            self.Nion = input['Nion']
        except Exception as e:
            self.logger.error("input dict must include Nion")
            self.logger.exception(e)
            raise
           
        try: 
            self.mi_amu = scipy.array(input['mi'])
            self.mi=self.mi_amu*scipy.constants.atomic_mass
        except Exception as e:
            self.logger.error("input dict must include mi")
            self.logger.exception(e)
            raise
                                
        try: 
            self.Params['f0'] = input['f0']
        except Exception as e:
            self.logger.error("input dict must include f0")
            self.logger.exception(e)
            raise

        try: 
            self.Params['B'] = input['B']
        except:
            self.Params['B'] = 0.0
            
        try: 
            self.Params['te'] = input['te']
        except:
            self.Params['te'] = 0.0
            
        try: 
            self.Params['alpha'] = input['alpha']
        except:
            self.Params['alpha'] = 0.0
            
        try: 
            self.Params['ne'] = input['ne']
        except:
            self.Params['ne'] = 0.0
            
        try: 
            self.Params['ti'] = scipy.array(input['ti'])
        except:
            self.Params['ti'] = scipy.zeros(self.mi.shape)
            
        try: 
            self.Params['ni'] = scipy.array(input['ni'])
        except:
            self.Params['ni'] = scipy.zeros(self.mi.shape)
            self.Params['ni'][0] = 1.0

        try: 
            self.Params['ve'] = input['ve']
        except:
            self.Params['ve'] = 0.0

        try: 
            self.Params['vi'] = scipy.array(input['vi'])
        except:
            self.Params['vi'] = scipy.zeros(self.mi.shape)

        try: 
            self.Params['nuen'] = input['nuen']            
        except:
            self.Params['nuen'] = 0.0

        try: 
            self.Params['nuin'] = scipy.array(input['nuin'])
        except:
            self.Params['nuin'] = scipy.zeros(self.mi.shape)                   

        if not self.override:
            self.checkParams()
        else:
            self.logger.info('Override flag is set, so not checking')
        
        return

    def checkParams(self):

        """
        checkParams: Checks parameters for shape consistency and composition
        
        """
        
        try:
            
            if not (
                    self.Params['ti'].shape[0] == self.Nion and 
                    self.Params['ni'].shape[0] == self.Nion and
                    self.Params['vi'].shape[0] == self.Nion and 
                    self.Params['nuin'].shape[0] == self.Nion and
                    self.mi.shape[0] == self.Nion 
                    ):
                raise ValueError('mui, ni, mi, vi, nuin must have length Nion')

            if scipy.absolute(1.0 - scipy.sum(self.Params['ni']))/len(self.Params['ni'])>1.0e-6:
                raise ValueError("ni must add to 1.0, currently adds to %2.7f" % scipy.sum(self.Params['ni']))
        
            for ion in self.mi_amu:
                if not str(int(ion)) in self.Amu2Ion.keys():
                    raise ValueError("ion masses must be one of %s" %self.Amu2Ion.keys())
                
        except Exception as e:
            self.logger.exception(e)
            raise
            
        return
        
    def computeSpec(self):

        """
        computeSpec: Does the job of computing the spectrum
        
        """
        
        start = time.time()
        #### compute parameters
        
        # magnetic aspect angle
        alpha = self.Params['alpha']*pi/180.0
        
        # if no Bfield, set alpha to 0 for ions and/or electrons
        if not self.emode[0]: # electrons
            alphae = 0.0
        else:
            alphae = alpha
        if not self.imode[0]: # ions
            alphai = 0.0
        else:
            alphai = alpha
        
        # k vector and wavelength (backscattering)
        k=2.0*pi/scipy.constants.c*(2.0*self.Params['f0']) 
        lambd=2.0*pi/k 
                    
        # electron and ion thermal speeds
        Ce=scipy.sqrt(scipy.constants.k*self.Params['te']/scipy.constants.m_e) 
        Ci=scipy.sqrt(scipy.constants.k*self.Params['ti']/self.mi)
        
        # electron and ion gyro frequencies
        Oe=scipy.constants.e*self.Params['B']/scipy.constants.m_e 
        Oi=scipy.constants.e*self.Params['B']/self.mi 
                
        # electron and ion plasma frequency
        wpe=scipy.sqrt(self.Params['ne']*scipy.constants.e**2/scipy.constants.m_e/scipy.constants.epsilon_0) 
        wpi=scipy.sqrt(self.Params['ni']*self.Params['ne']*scipy.constants.e**2/self.mi/scipy.constants.epsilon_0)
        
        # electron, ion, and plasma debye length
        he=Ce/wpe 
        
        # coulomb collisions
        if self.emode[1]:
            nuei = 1e-6*54.5*self.Params['ne']/self.Params['te']**1.5 # Schunk and Nagy 4.144
            nuee = 1e-6*54.5/scipy.sqrt(2.0)*self.Params['ne']/self.Params['te']**1.5 # Schunk and Nagy 4.145
            
            psiec_par = nuei/(k*Ce*scipy.sqrt(2.0))
            psiec_perp = (nuei+nuee)/(k*Ce*scipy.sqrt(2.0))
                                               
        else:
            psiec_par = 0.0
            psiec_perp = 0.0
            
        if self.imode[1]:
            psiic = scipy.zeros(self.mi.shape)        
            for i1 in range(self.Nion):
                for i2 in range(self.Nion):
                    if self.Params['ni'][i2]>1.0e-6:
                        tb = collisionCoeffs.getBst(self.mi_amu[i1],self.mi_amu[i2])
                        psiic[i1] = psiic[i1] + tb*self.Params['ni'][i2]*self.Params['ne']*1.0e-6/self.Params['ti'][i2]**1.5 # Schunk and Nagy 4.143                        
            psiic = psiic/(k*Ci*scipy.sqrt(2.0))
        else:
            psiic = scipy.zeros(self.mi.shape)
        
        # neutral collisions
        if self.emode[2]:
            psien = self.Params['nuen']/(k*Ce*scipy.sqrt(2))
        else:
            psien = 0.0
        if self.imode[2]:
            psiin = self.Params['nuin']/(k*Ci*scipy.sqrt(2))
        else:
            psiin = scipy.zeros(self.mi.shape)
                
        ####
 
        # frequency        
        if self.fmax>0.0:
            fmax = self.fmax
        else:
            fmax=scipy.sum(self.Params['ni']*k*Ci*scipy.sqrt(2.0))/2.0+scipy.absolute(self.Params['vi']).max()/lambd
        
        fmin=-fmax
        ff=scipy.linspace(fmin,fmax,self.Nfreq);
        ww=2.0*pi*ff
        
        # special normalized frequencies
        thetae = (ww + k*self.Params['ve'])/(k*Ce*scipy.sqrt(2.0)) 
        thetai=(scipy.repeat(ww[:,scipy.newaxis],self.Nion,axis=1) + k*self.Params['vi'])/(k*Ci*scipy.sqrt(2.0))  
        
        phie=Oe/(k*Ce*scipy.sqrt(2))
        phii=Oi/(k*Ci*scipy.sqrt(2))
        
        # compute Gordeyev integrals
        Ji=scipy.zeros(thetai.shape,dtype='complex64')
        Je=scipy.zeros(thetae.shape,dtype='complex64')    
        
        Je,Nl=self.gordeyev_cz(thetae,(alphae,phie,psien,psiec_par,psiec_perp),czparams=self.czparamse)

        #print Nl
        for ithi in range(self.Nion):
            if self.Params['ni'][ithi]>1.0e-6:
                Ji[:,ithi],Nl=self.gordeyev_cz(thetai[:,ithi],(alphai,phii[ithi],psiin[ithi],psiic[ithi],psiic[ithi]),czparams=self.czparamsi)        
                
        # compute spectrum
        spec=self.compute_spec(Je,Ji,thetae,thetai,Ce,Ci,self.Params['ni'],self.Params['te']/self.Params['ti'],he,k)
        spec = spec*self.Params['ne']        

        # compute the acf
        tau,acf=self.spec2acf(ff,spec)
       
        # normalize if requested
        if self.normalize==1:
            acf=acf.real
            acf=acf/acf[int(scipy.floor(tau.size/2.0))]
    
        return ff,spec,tau,acf
            
    def compute_spec(self,Je,Ji,thetae,thetai,Ce,Ci,ni,mui,he,k):
        
        # ions
        sis = scipy.zeros(thetae.shape,dtype='complex64')
        nthi = scipy.zeros(thetae.shape,dtype='float64')
        for ithi in range(self.Nion):
            # ion conductivity
            sis = sis + ni[ithi]*mui[ithi]*(1.0 - 1.0j*thetai[:,ithi]*Ji[:,ithi])/(k*k*he*he)
            # ion density thermal fluctuation spectrum
            nthi = nthi + ni[ithi]/scipy.sqrt(2.0)/k/Ci[ithi]*2.0*Ji[:,ithi].real
            
        # electron conductivity
        ses = (1.0 - 1.0j*thetae*Je)/(k*k*he*he)

        # electron thermal fluctuation spectrum
        nthe = 2.0*Je.real/scipy.sqrt(2.0)/k/Ce
        
        # numerators and denominator
        num1 = scipy.power(scipy.absolute(1.0+sis),2.0)*nthe
        num2 = scipy.power(scipy.absolute(ses),2.0)*nthi
        den = scipy.power(scipy.absolute(1.0+sis+ses),2.0)

        # spectrum
        spec=(num1+num2)/den
                
        return spec	        
               
               
    # gordeyev
    def gordeyev_cz(self,theta,gvparams,czparams=(1e-6,2e5,100,10,1.0)):
        '''
            Computes Gordeyev integral for either electrons or ions
            Based on Chirp-z transform method of Li et al., IEEE Trans. Ant. Prop., 1991
        
            Inputs:
                theta - normalized frequency, (w - kV)/(k*Vth)
                gvparams - list containing:
                    alpha - magnetic field aspect angle, radians
                    phi - normalized gyro frquency, O/(k*Vth)
                    psin - normalized neutral collision frquency, nu/(k*Vth)
                    psic_par - normalized coulomb collision frquency parallel, nu/(k*Vth)
                    psic_perp - normalized coulomb collision frquency perpendicular, nu/(k*Vth)
            Optional Inputs:
                czparams - list containing:
                    tol - tolerance for integral
                    Nmax - maximum number of points in integral
                    maxLoops - maximum number of integral loops
                    Nstart - starting number of points in integral
                    kmax - upper bound of integral (starting value)
        '''
        
        (alpha,phi,psin,psic_par,psic_perp) = gvparams
        (tol,Nmax,maxLoops,Nstart,kmax) = czparams
        Sn, flag, loop = sommerfeldIntegralM(self.funcgv,theta,0.0,bs=kmax,Ns=Nstart,Nmax=Nmax,maxLoops=maxLoops,tol=tol,additParams=gvparams)
        
        #print loop
        J = Sn/(1.0-psin*Sn) 

        return J,loop

    def funcgv(self,t,inputParams):

        """
        funcgv: returns Gordeyev integrand for case of magnetic field, coulomb collisions, neutral collisions

        input parameters are a list containing:
            alpha - magnetic field aspect angle, radians
            phi - normalized gyro frquency, O/(k*Vth)
            psin - normalized neutral collision frquency, nu/(k*Vth)
            psic_par - normalized coulomb collision frquency parallel, nu/(k*Vth)
            psic_perp - normalized coulomb collision frquency perpendicular, nu/(k*Vth)
        
        """
        
        (alpha,phi,psin,psic_par,psic_perp) = inputParams
        
        # collisional term
        K0 = scipy.exp(-t*psin) 
        
        # coulomb collisions
        if (psic_par+psic_perp)>1e-6:
            
            alpha=pi/2.0-alpha

            ca2 = scipy.power(scipy.cos(alpha),2.0)
            sa2 = scipy.power(scipy.sin(alpha),2.0)

            # parallel term
            K1 = scipy.exp(-(0.5*sa2/(psic_par*psic_par)*(psic_par*t-1.0+scipy.exp(-psic_par*t)))) 

            # perpendicular term
            gamma = scipy.arctan(psic_perp/phi)
            K2 = scipy.exp(-(0.5*ca2/(phi*phi+psic_perp*psic_perp)*(scipy.cos(2.0*gamma)+psic_perp*t-scipy.exp(-psic_perp*t)*scipy.cos(phi*t-2.0*gamma)))) 

        else:

            ca2 = scipy.power(scipy.cos(alpha),2.0)
            sa2 = scipy.power(scipy.sin(alpha),2.0)
            
            K1=scipy.exp(-sa2*scipy.power(scipy.sin(0.5*phi*t),2.0)/(phi*phi)-0.25*t*t*ca2)
            K2=1.0        

        val = K0*K1*K2
        

        return val  
        
    def spec2acf(self,f,s):
        # converts the spectra to an acf

        Nspec=s.size

        zsize=Nspec*3
        if scipy.mod(Nspec,2.0)==0.0:
            zsize=Nspec/2.0

        spec=scipy.concatenate((scipy.zeros((zsize),dtype='float64'),s,scipy.zeros((zsize),dtype='float64')),axis=0) # zero pad the spectra
        
        NFFT=spec.size
        df=f[1]-f[0]
        tau=scipy.linspace(-1.0/df/2.0,1.0/df/2.0,NFFT)
        dtau=tau[1]-tau[0]

        m=scipy.fftpack.fftshift(scipy.fftpack.ifft(scipy.fftpack.ifftshift(spec)))/dtau
                
        return tau, m
    