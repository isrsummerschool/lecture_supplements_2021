"""
CoordConv.py

Written by: Patrick Sullivan, Jarred Velazauez, Landon Airey
Updated by: Michael Nicolls

This is a Coordinate Conversion library of helper functions
for the AMISR.

"""

import sys
import numpy as np
import scipy
from scipy import linalg

import numpy as np
pi = np.pi

def elaz2Dir(elazPair, degFlag = False):
    """
    ##########################################################
    # Method:elaz2Dir #                                      #
    ###################                                      #
    # Inputs:   elazPair, degFlag                            #
    # Outputs:  dirPair                                      #
    #                                                        #
    # Summary:  The method will convert azimuth-elevation    #
    #           coordinates to directional sine coordinates  #
    #           It returns a list of corresponding thx-thy   #
    #           coordinates                                  #
    #           It expects a list of tuples entered as:      #
    #           listOfelaz = [(el0, az0),(el1,az1)...,       #
    #                        (eln,azn)]                      #
    #           If degFlag is equal to 1, the input          #
    #           angles are in degrees.                       #
    #                                                        #
    ##########################################################
    """      
   
    dirPair = []
            
    for i in range(len(elazPair)):
        (el,az) = elazPair[i]

        if degFlag == True:
            el = el*pi/180
            az = az*pi/180

        x = np.cos(el)*np.cos(az)
        y = np.cos(el)*np.sin(az)
        z = np.sqrt(1.0-(x**2.0+y**2.0))
        
        thx = np.around(np.arcsin(x), decimals = 10)
        thy = np.around(np.arcsin(y), decimals = 10)
        
        if degFlag == True:
            thx = thx*180.0/pi
            thy = thy*180.0/pi
        
        dirPair.append((thx,thy))
    
    return dirPair
    
def dir2elaz(dirPair, degFlag = False, tilt = (0.0,0.0)):
    """
    ##########################################################
    # Method:dir2elaz #                                      #
    ###################                                      #
    # Inputs:   dirPair, degFlag                             #
    # Outputs:  elazPair (in radians)                        #
    #                                                        #
    # Summary:  The method will convert directional sine     #
    #           coordinates to azimuth-elevation coordinates.#
    #           It returns a list of corresponding el-az     #
    #           coordinates                                  #
    #           It expects a list of tuples entered as:      #
    #           listOfdir = [(thx0, thy0),(thx1,thy1)...,    #
    #                        (thxn,thyn)]                    #
    #           If degFlag is equal to 1, the input          #
    #           angles are in degrees.                       #
    #                                                        #
    ##########################################################
    """

    R_GF, R_FG = Rotate_Face(tilt)
    elazPair = []
 
    for i in range(len(dirPair)):
        (thx,thy) = dirPair[i]
        
        if degFlag == True:
            xhat = np.sin(thx*pi/180)
            yhat = np.sin(thy*pi/180)
            zhat = np.sqrt(1.0-xhat**2.0-yhat**2.0)    
        
        else:
            xhat = np.sin(thx)
            yhat = np.sin(thy)
            zhat = np.sqrt(1.0-xhat**2.0-yhat**2.0)

        x = xhat*R_FG[0,0] + yhat*R_FG[0,1] + zhat*R_FG[0,2]
        y = xhat*R_FG[1,0] + yhat*R_FG[1,1] + zhat*R_FG[1,2]
        z = xhat*R_FG[2,0] + yhat*R_FG[2,1] + zhat*R_FG[2,2]

        el = np.arcsin(z)
        az = np.arctan2(-y,x)

        if degFlag == True:
            el = el*180.0/pi
            az = az*180.0/pi

        elazPair.append((el,az))        

    return elazPair

def geo2faceElAz(geoElAzPair, degFlag = False, tilt = (0.0,0.0)):
    """
    ##########################################################
    # Method:geo2faceElAz #                                  #
    #######################                                  #
    # Inputs:   geoElAzPair (tuple list of geodetic El and Az#
    #               coordinates [(El,Az), ... ,(El, Az)] )   #
    #           degFlag (Set to True if ElAz are in degrees) #
    #                                                        #
    # Outputs:  faceElAzPair (tuple list of local El and Az  #
    #               coordinates [(El,Az), ... ,(El, Az)] )   #
    #                                                        #
    # Summary:  Takes in geodetic relative Azimuth Elevation #
    #           coordinates and converts them to face local  #
    #           Azimuth Elevation coordinates                #
    #                                                        #
    ##########################################################
    """
    R_GF, R_FG = Rotate_Face(tilt)
    faceElAzPair = []
 
    for i in range(len(geoElAzPair)):
        (geoEl, geoAz) = geoElAzPair[i]
        
        if degFlag == True:
            geoEl = geoEl*pi/180.
            geoAz = geoAz*pi/180.

        #Define xyz coords (North, West, Up) for geodetic coords
        xg = np.cos(geoEl)*np.cos(geoAz)
        yg = -np.cos(geoEl)*np.sin(geoAz)
        zg = np.sqrt(1.0 - (xg**2.0 + yg**2.0))
        
        #Convert geodetic to face relative xyz coords using the 
        # R_GF transformation matrix
        xf = xg*R_GF[0,0] + yg*R_GF[0,1] + zg*R_GF[0,2]
        yf = xg*R_GF[1,0] + yg*R_GF[1,1] + zg*R_GF[1,2]
        zf = xg*R_GF[2,0] + yg*R_GF[2,1] + zg*R_GF[2,2]

        faceEl = scipy.arcsin(zf)
        faceAz =  scipy.arctan2(yf,xf)
        
        if degFlag == True:
            faceEl = faceEl*180.0/pi
            faceAz = faceAz*180.0/pi

        faceElAzPair.append((faceEl,faceAz))     

    return faceElAzPair


def Rotate_Face(tilt):
    """
    ##########################################################
    # Rotate_Face  #                                         #
    ################                                         #
    # Inputs:   N/A                                          #
    #                                                        #
    # Outputs:  R_GF (Trans. matrix from Geodetic to Face)   #
    #           R_FG (Trans. matrix from Face to Geodetic)   #
    #                                                        #
    # Summary:  This helper function defines a transformation#
    #           matrix to convert x,y,z coordinates in face  #
    #           relative coordinates to geodetic relative    #
    #           x,y,z coordinates, or vice versa.            #
    #                                                        #
    ##########################################################
    """
    AZ_ROT = scipy.deg2rad(tilt[1])
    E_TILT = scipy.deg2rad(tilt[0])

    #R_FG1 describes step 1 (transform xyz geodetic to a new xyz)
    R_FG1 = scipy.array([[ scipy.cos(AZ_ROT),scipy.sin(AZ_ROT), 0 ],
                         [ -scipy.sin(AZ_ROT), scipy.cos(AZ_ROT), 0 ],
                         [ 0, 0, 1 ]])
    #R_FG2 describes step 2 (transform new xyz to xyz face relative)
    R_FG2 = scipy.array([[ scipy.cos(E_TILT), 0, scipy.sin(E_TILT)],
                         [ 0, 1, 0 ],
                         [ -scipy.sin(E_TILT), 0, scipy.cos(E_TILT)]])

    #R_FG is the total transformation
    R_FG = R_FG1.dot(R_FG2)

    #R_FG is the inverse of the R_GF transformation and 
    # defines a face to geodetic transformation    
    R_GF = scipy.linalg.inv(R_FG) 

    return R_GF, R_FG
