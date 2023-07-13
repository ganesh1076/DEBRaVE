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
import argparse
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits


#-------------------------------global vars-----------------------------------#

# Verbose mode placeholder
verboseprint = lambda *a, **k: None

#--------------------------------classes--------------------------------------#

class Spectra:
    """
    This is the superclass for all spectra objects.

    Attributes:
        time (str): the timestamp of the observation.
        radec (tuple): the RA and DEC of the observed object.
        spectra_data (numpy array): the wavelength and flux data of the spectra.
            2D array of shape (N,2).
        rms (float): the root mean square of the measured fluxes.
        errors (list): a list of error messages.

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

        Args:
            error_msg (str): the error message to be added to the error list.
        
        Returns:
            None

        """

        self.errors.append(error_msg)
        print(error_msg)


    def plotSpectra(self, savename=None):
        """
        Plots the spectra. Saves the plot if a savename is provided.
        Returns the matplotlib figure object.

        Args:
            savename (str, optional): the name of the file to save the plot to if provided.
                Defaults to None.
        
        Returns:
            matplotlib figure: the figure object of the plot.

        """

        # Plot the spectra
        plt.plot(self.spectra_data[:,0], self.spectra_data[:,1])
        plt.xscale('log')
        plt.gca().invert_xaxis()

        # Adjust the plot parameters
        plt.title(f"Spectra for {self.time}\nRA: {self.radec[0]:.1f}, DEC: {self.radec[1]:.1f}")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Flux (arbitrary units)")
        plt.grid(True)

        # Save the plot if a savename is provided
        if (savename is not None):
            plt.savefig(savename)

        # Return the figure object
        return plt.gcf()


class FStarSpectra(Spectra):
    """
    This is the class for F-star spectra.
    """

    def __init__(self, time, radec, wavelength, fluxes):
        super().__init__(time, radec, wavelength, fluxes)

    
    def crossCorrelate(self, template):
        """
        Cross-correlates the spectra with a template spectra.
        *Note: the wavelengths in the two spectra must be the same.

        Args:
            template (FStarSpectra): the template spectra to cross-correlate with.

        Returns:
            numpy array: the cross-correlation function.

        """

        if len(self.spectra_data) != len(template.spectra_data):
            self.addError("Primary and secondary spectra are of different lengths!")
            return
        elif (self.spectra_data[:,0] != template.spectra_data[:,0]).all():
            self.addError("Primary and secondary spectra have different wavelengths!")
            return

        # Define the cross correlation variables
        N = len(self.spectra_data)
        sig_g = self.rms
        sig_t = template.rms

        # Calculate the cross correlation
        cross_corr = (N*sig_g*sig_t)**-1 * np.convolve(self.spectra_data[:,1], template.spectra_data[:,1])
        return cross_corr


    def mapTODCOR(self, template1, template2, light_ratio=1, savename=None, plot_block=True):
        """
        Make a heatmap of the TODCOR cross correlation function.

        Args:
            template1 (FStarSpectra): the primary template spectra.
            template2 (FStarSpectra): the secondary template spectra.
            light_ratio (float, optional): the ratio of the light from the primary to the secondary.
        
        Returns:
            numpy array: the TODCOR index array.

        """

        # Obtain TODCOR index array
        ind_arr = self.TODCOR(template1, template2, light_ratio)

        # Plot the heatmap
        plt.pcolormesh([template1.spectra_data[:,0], template2.spectra_data[:,0]], ind_arr, cmap='seismic')
        plt.colorbar()
        plt.title("TODCOR Cross-Correlation Function")
        plt.xlabel("Unshifted Wavelength (nm)")
        plt.ylabel("Shifted Wavelength (nm)")

        # Save the heatmap if a savename is provided
        if (savename is not None):
            plt.savefig(savename)

        plt.show(block=plot_block)
        return ind_arr
    

    def TODCOR(self, template1, template2, light_ratio=1):
        """
        Performs the TODCOR algorithm on the spectra.
        
        Args:
            template1 (FStarSpectra): the primary template spectra.
            template2 (FStarSpectra): the secondary template spectra.
            light_ratio (float, optional): the ratio of the light from the primary to the secondary.
                Defaults to 1.
        
        Returns:
            numpy array: the TODCOR index array.
        """

        # Obtain individual cross correlation functions
        cross_corr1 = self.crossCorrelate(template1)
        cross_corr2 = self.crossCorrelate(template2)
        cross_corr12 = template1.crossCorrelate(template2)

        # Calculate the TODCOR cross correlation function
        cross_corr = np.empty((len(template1.spectra_data),len(template2.spectra_data)))
        for i in range(len(template1.spectra_data)):
            for j in range(len(template2.spectra_data)):
                cross_corr[i,j] = (cross_corr1[i] + light_ratio*cross_corr2[j])/np.sqrt(1 + 2*light_ratio*cross_corr12[i,j] + light_ratio**2)


        return cross_corr


#-------------------------------functions-------------------------------------#

def readSpectraFITS(filename):
    """
    Reads the header and contents of a spectra file. Returns the header
    as a dictionary and the body as a numpy object.

    Args:
        filename (str): the name of the file to be read.
    
    Returns:
        dict: the header of the file as a dictionary.
        numpy.ndarray: the body of the file as a numpy object.

    """
    verboseprint(f"Unpacking {filename}...")

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

    verboseprint(f"Header fields: {list(header.keys())}")

    return time, (ra, dec), wavelengths, fluxes



#---------------------------------main----------------------------------------#


def main(fits_files, templates, light_ratio=1, parallel=False):

    # Initialize the template spectra objects
    # TODO: make this more specific to the template formats
    verboseprint(f"{templates[0]}\n{templates[1]}")
    template1 = FStarSpectra(*readSpectraFITS(templates[0]))
    template2 = FStarSpectra(*readSpectraFITS(templates[1]))

    # Initialize the spectra objects and evaluate the cross correlation
    # for each. Then, plot the cross correlation.
    if (parallel):
        # Obtain the number of cores for multiprocessing
        num_cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(num_cores)

        # Define arguments extraction and mapping functions from FITS to Specta objects
        def processFITS(filename):
            star_spectrum = FStarSpectra(*readSpectraFITS(filename))
            return star_spectrum.mapTODCOR(template1, template2, light_ratio, plot_block=False)

        # Map the arguments to the multiprocessing pool
        try:
            TODCOR_ind_maps = pool.starmap(processFITS, fits_files)
        except Exception as err:
            print(f"Error: {err}")

        # Close the multiprocessing pool
        pool.close()
        pool.join()


    else:
        # Iterate through the spectrum files and return the TODCOR cross correlation
        for spectrum_file in fits_files:
            star_spectrum = FStarSpectra(*readSpectraFITS(spectrum_file))
            star_spectrum.mapTODCOR(template1, template2, light_ratio)

    pass


if __name__ == "__main__":

    # Generate argument parser
    arg_parser = argparse.ArgumentParser(description="DEBRaVE: Double Eclipsing Binary Research and Visualization Engine",
                                         formatter_class=argparse.RawTextHelpFormatter)

    # Add arguments
    arg_parser.add_argument("-s", "--spectra", type=str, nargs='+', required=True,
                            help="F-Class star spectra to be read from FITS file(s).")
    arg_parser.add_argument("-t", "--templates", type=str, nargs=2, required=True,
                            help="Template FITS files to be read. Two required.")
    arg_parser.add_argument("-l", "--light_ratio", type=float, default=1,
                            help="Scaling light ratio of the unshifted to shifted templates.")
    arg_parser.add_argument("-p", "--parallel", action="store_true",
                            help="Run the program in parallel mode.")
    arg_parser.add_argument("-v", "--verbose", action="store_true",
                            help="Run the program in verbose mode.")

    # Parse arguments
    args = arg_parser.parse_args() # cml argument dictionary
    fits_files = args.spectra
    templates = args.templates
    light_ratio = args.light_ratio
    parallel = args.parallel

    # Check for verbose mode
    if (args.verbose):
        verboseprint = print

    main(fits_files, templates, light_ratio, parallel)