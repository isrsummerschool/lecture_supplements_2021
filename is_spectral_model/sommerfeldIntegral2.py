
"""

Routines to evaluate Sommerfeld integral (i.e., Gordeyev integrals)

M Nicolls - 2013

"""

import logging
import scipy, scipy.special,numpy
import time

pi = scipy.pi

def sommerfeldIntegralM(funct,x,a,bs=1.0,Ns=10,Nmax=1e4,Nmin=10,tol=1e-6,maxLoops=100,additParams=()):

    '''
        Computes Sommerfeld-type integral of the form:
            Integral from a to b of f(k)*exp(alph*k*x) dk
            
        Based on Erf transform method of Ooi, B.L., Microw. Opt. Tech. Lett., 2007
    
        Inputs:
            funct - function handle that evaluates integrand and takes two inputs: x and a list of optional parameters
            x - range values to evaluate integral at
            a - starting value of integral
        Optional Inputs:
            bs - ending value of integral (will be adjusted as it iterates)
            Ns - starting number of points in integral (actually will be 2xNs+1)
            Nmax - maximum number of points allowed
            Nmin - minimum number of points allowed
            tol - desired tolerance
            maxLoops - maximum number of loops (iterations) allowed
            additParams - list of additional parameters to pass to funct
            
        Outputs:
            (S, flag, loop)
            S - integral at points defined in x
            flag - 1 or -1 if successful or failed, respectively
            loop - loop index at termination
            
        Note - future improvements: 
                - could speed up using previous value / Simpson's rule
                - need better way of scaling and choosing bs - can be slow to converge
                - likely becomes much slower that chirp-z method for large number of x values
    '''
    
    logger = logging.getLogger(__name__)
    
    # fixed variables used in algorithm
    initAssymFuncVal = 0.05
    initFuncDiffVal = 0.2
    boundsScalar = 1.0
    npointsScalar = 1.5
    
    #
    Nx = x.shape[0]
    N=max([Ns,Nmin])
    b=bs

    #
    Sp = 1e10
    loop = 0
    while 1:

        # abort if N is getting too large
        if N>Nmax or loop>maxLoops:
            logger.error('Max N or Max loops reached; aborting')
            flag = -1
            break   
                
        #
        h = 1.0/N*scipy.log(1.05*scipy.sqrt(2.0*N))    
        n = scipy.arange(-N,N+1)

        #
        An = scipy.cosh(n*h)*scipy.exp(-scipy.power(scipy.sinh(n*h),2.0))
        inval = 0.5*(b+a)+0.5*(b-a)*scipy.special.erf(scipy.sinh(n*h))

        # integrand evaluation
        f = funct(inval,additParams)

        # if first time through, we do a check to make sure integration bounds
        # and number of points are reasonable
        if loop == 0:
            # refine end value of integral 
            # note this is adjusted in each iteration by boundsScalar 

            if scipy.mean(f[-Nmin:])>initAssymFuncVal:
                b*=2
                continue
            
            # after that, refine number of points to capture variation
            # note this is adjusted in each iteration by npointsScalar 
            if scipy.absolute(scipy.diff(f)).max()>initFuncDiffVal:
                N=N*2
                continue

        # array Nx rows by 1 col
        Amat = An*f

        # matrix Nx rows by N cols
        
        x1 = scipy.dot(inval[:,scipy.newaxis],x[scipy.newaxis,:])

        fkern = scipy.exp(-1.0j*x1)

        #
        S = h*(b-a)/scipy.sqrt(pi)*scipy.dot(Amat,fkern)
        
        # compute mean squared error
        ds = scipy.sum(scipy.power(scipy.absolute((S - Sp)/Sp),2.0))/Nx
        

        # compute maximum absolute error
        # ds = scipy.absolute((S - Sp))/scipy.absolute(Sp)
        # ds = ds.max()
        
        if ds<tol:
            flag=1
            break

        # increment b and number of points
        Sp = S[:]
        N *= npointsScalar
        b *= boundsScalar
        loop += 1

    return S, flag, loop

def sommerfeldIntegral(funct,x,a,bs=1.0,Ns=10,Nmax=1e4,tol=1e-6,maxLoops=100,additParams=()):
    
    '''
        Same as sommerfeldIntegralM but does not handle array inputs (x must be a scalar)
    '''
    
    logging.getLogger(__name__)
    
    N=Ns
    b=bs
    
    Sp = 1e10
    while 1:
                
        h = 1.0/N*scipy.log(1.05*scipy.sqrt(2.0)*N)    
        n = scipy.arange(-N,N+1)
        
        An = scipy.cosh(n*h)*scipy.exp(-scipy.power(scipy.sinh(n*h),2.0))
        inval = 0.5*(b+a)+0.5*(b-a)*scipy.special.erf(scipy.sinh(n*h))
        
        f=funct(inval,additParams)*scipy.exp(-1.0j*x*inval)

        S = h*(b-a)/scipy.sqrt(pi)*scipy.sum(An*f)

        ds = scipy.absolute(S-Sp)/scipy.absolute(Sp)
                
        if ds<tol:
            flag=1
            break

        # abort if N is getting too large
        if N>Nmax or loop>maxLoops:
            logger.error('Max N or Max loops reached; aborting')
            flag = -1
            break    
            
        Sp = S*1.0
        N = N*2
        b = b*1.5

    return S, flag

def funcEg(k,a):
        
    return scipy.exp(-k**2.0/2)

if __name__ == "__main__":
    '''
    This example taken from Drachman et al., IEEE, 1989
    '''

    # test 1
    
    x = scipy.arange(2)                
    out, flag, loop = sommerfeldIntegralM(funcEg,x,0.0,bs=1.0,Ns=10,Nmax=1e4,tol=1e-6)
    print(out[1].real*2.0/scipy.sqrt(2.0*pi) - scipy.exp(-0.5))
    
    print(out)
