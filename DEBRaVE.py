"""
Filename:   DEBRaVE.py
Author(s):  Peter Quigley, David Dougan, Ganesh Pawar
Contact:    pquigley@uwo.ca
Created:    2023-07-10
Updated:    2023-07-10
    
Usage: python DEBRaVE.py
"""

# Module imports
import os, sys
import astropy
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits

#-------------------------------global vars-----------------------------------#


#--------------------------------classes--------------------------------------#

class Spectra:
    """
    This is the superclass for all spectra objects.
    """

    def __init__(self, time, radec, wavelength, fluxes):

        # Error collection
        self.errors = []

        # Parse init arguments
        self.time = time  # timestamp of the spectra as a string
        self.radec = radec  # this should be a tuple of (ra, dec)
        self.spectra_data = np.array([wavelength, fluxes]).transpose()
        self.rms = np.sqrt(np.mean(fluxes**2))  # root mean square of the fluxes

    def addError(self, error_msg):
        """
        Error handling function for spectra objects.
        """

        self.errors.append(error_msg)
        print(error_msg)


class FStarSpectra(Spectra):
    """
    This is the class for F-star spectra.
    """

    def __init__(self, time, spectra_data):
        super().__init__(time, spectra_data)

    
    def crossCorrelate(self, template):
        """
        Cross-correlates the spectra with a template spectra.
        """

        if len(self.spectra_data) != len(template.spectra_data):
            self.addError("Real and refence spectra are of different lengths!")

        # Define the cross correlation variables
        N = len(self.spectra_data)
        sig_g = self.rms
        sig_t = template.rms

        # Calculate the cross correlation
        cross_corr = (N*sig_g*sig_t)**-1 * np.convolve(self.spectra_data, template.spectra_data)


#-------------------------------functions-------------------------------------#

def readSpectraFITS(filename):
    """
    Reads the header and contents of a spectra file. Returns the header
    as a dictionary and the body as a numpy object.
    """

    with fits.open(filename) as hdul:
        # Read the header
        header = hdul[0].header
        time = header['DATE-OBS']
        ra   = header['RA']
        dec  = header['DEC']

        # Read the body
        body = hdul[1].data[0]
        wavelengths = body[0]
        fluxes = body[1]

    print(f"Header fields: {header.keys()}")

    return header, wavelengths, fluxes


#---------------------------------main----------------------------------------#


def main():

    pass