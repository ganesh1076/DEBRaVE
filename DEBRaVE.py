"""
Filename:   DEBRaVE.py
Author(s):  Peter Quigley, David Dougan, Ganesh Pawar
Contact:    pquigley@uwo.ca
Created:    2023-07-10
Updated:    2023-07-10
    
Usage: python DEBRaVE.py [-d][-p][-r][-s]
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

    def __init__(self):

        self.errors = []
        pass

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

    def __init__(self):
        pass

#-------------------------------functions-------------------------------------#


#---------------------------------main----------------------------------------#
