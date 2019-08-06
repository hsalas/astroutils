#! /usr/bin/python
import numpy as np
from astropy import units as u
from astropy.constants import c
from astropy.io import fits

# to have auto complete an navigate directories ofr the user input
import readline
readline.parse_and_bind('tab: complete')
readline.set_completer_delims(' \t\n')

# Functions for do the unit conversion


def cfccd_to_mjy_per_pix(image):
    """Convert from CTIO 0.9m telescope taken with cfccd camera to mJy/px

    inputs
            image: image to be converted
            band: (str) FUV or NUV
    """
    # zero point magnitudes
    f_zp = {'U':1884, 'B':4646, 'V':3953, 'R':2875}
    # read values from headers
    band = image[0].header['FILTER']
    zp = image[0].header['ZP']

    # do conversion
    data_mjy = image[0].data * f_zp[band] * 10**(-0.4*zp) * 10**3
    image[0].data = data_mjy

    # update header
    image[0].header['UNITS'] = 'mJy/pixel'
    image[0].header['history'] = 'units converted to  mJy/pixel'

    return image


def galex_to_mjy_per_pix_2(image, band):
    """Convert GALEX images (FUV or NUV) from CPS (counts per second) per pixel
    to mJy per pixel. Using the zero point values for the GALEX bands

    Inputs
            image: Galex image in CPS/pix
            band: (str) FUV or NUV
    Output
            image: Galex image in mJy/pix

    Formulas
            F[Jy] = f*3631*10**(-0.4*zp)
    """

    band = band.upper()

    # instrument constants
    zp = {"FUV": 18.82, "NUV": 20.08}
    # convert from CPS to flux (mJy)
    image[0].data *= 3631 * 10**(-.4*zp[band]) * 1e3
    # update header
    image[0].header['UNIT'] = 'mJy/pixel'
    image[0].header['history'] = 'units converted to  mJy/pixel'

    return image


def galex_to_mjy_per_pix(image, band):
    """Convert GALEX images (FUV or NUV) from CPS (counts per second) per pixel
    to mJy per pixel. The flux is converted to F_lambda using the formula given
    below and then is F_nu and into mJy.

    Inputs
            image: Galex image in CPS/pix
            band: (str) FUV or NUV
    Output
            image: Galex image in mJy/pix

    Formulas
            FUV: Flux [erg sec-1 cm-2 Angstrom-1] = 1.40 x 10-15 x CPS
            NUV: Flux [erg sec-1 cm-2 Angstrom-1] = 2.06 x 10-16 x CPS
    """

    band = band.upper()

    # instrument constants
    conv_factor = {"FUV": 1.40e-15, "NUV": 2.06e-16}
    eff_wl = {"FUV": 1528.0 * u.AA, "NUV": 2271.0 * u.AA}

    # convert to f_lambda
    image_flambda = image[0].data * conv_factor[band] * (u.erg / u.s / u.cm **
                                                         2 / u.AA)

    # convert to f_nu
    image_fnu = (image_flambda * (eff_wl[band]**2 / c)).to(u.erg / u.s / u.cm
                                                           ** 2 / u.Hz)

    # convert to mJy
    image_mJy = image_fnu.to(u.mJy)
    image[0].data = image_mJy.value

    # update header
    image[0].header['BUNIT'] = 'mJy/pixel'
    image[0].header['history'] = 'units converted to  mJy/pixel'

    return image


def f_2mass_to_mj_per_pix(image, band):
    """Convert 2MASS images (J, H or Ks) in DN (dta number units) to mJy/px.
    The converstion is perform using the zeropoint valus and zero point
    magnitudes for each of the 2MASS bands.

    Inputs
            image: 2MASS image
            band: (str) J, H or Ks
    Output
            image: 2MASS image in mJy/pix
    Formulas:

            i)  mag = m_zp -2.5 log10(x)
            ii) flux = flux_zp*10^(-mag/2.5)

            from i) and ii) ==>

            iii) flux = flux_zp*10^(-mag/2.5)*x
    """

    band = band.upper()

    # zero poitns
    flux_zp = {'J': 1594.0, 'H': 1024.0, 'K': 666.7}
    mag_zp = {'J': 21.02109909, 'H': 20.48749924, 'K': 19.96899986}

    # to flux
    image_Jy = flux_zp[band]*10**(-.4*mag_zp[band])*image[0].data*u.Jy
    image_mJy = image_Jy.to(u.mJy)
    image[0].data = image_mJy.value

    # update header
    image[0].header['BUNIT'] = 'mJy/pixel'
    image[0].header['history'] = 'units converted to  mJy/pixel'

    return image


def spitzer_to_mJy_per_pix(image, band):
    """
    Converts Spitzer photometric images (IRAC and MIPS) from MJy/sr to mJy/pix

    Inputs
        image: Spitzer image in MJy/sr units
        band: (str) irac1, irac2, irac3, irac4, mips24, mips70 or mips160
    Output:
        image: Spitzer converted to  mJy/pix units
    """

    # get value from header
    pix_x = np.abs(image[0].header['CD1_1'])
    pix_y = np.abs(image[0].header['CD2_2'])

    # get conv factor from str to pix
    conv_factor = pix_x * pix_y * 3600 * 3600 / u.sr.to(u.arcsec**2) * 1e9
    image_mJy = image[0].data * conv_factor
    image[0].data = image_mJy

    # irac aperture corretion
    # ap_corr = {'irac1':0.91, 'irac2':0.94, 'irac3':0.73, 'irac4':0.74}
    # if 'irac' in band:
    #    image[0].data *= ap_corr{band}

    # update header
    image[0].header['BUNIT'] = 'mJy/pixel'
    image[0].header['history'] = 'units converted to  mJy/pixel'
    return image


def herschel_to_mJy_per_pix(image, band):
    """Converts units from Herschel photometric images (PACS and SPIRE) from Jy
    per pixel to mJy per pixel.

        Inputs
            image: Herschel image in Jy/pix units

        Output:
            image: Herschel converted to  mJy/pix units
    """

    # from Jy to mJy
    image_mJy = image[0].data * 1e3

    # factor to update phometry
    # conv_factor = {70:1.05, 100:1.08, 160:1.0, 250:.92, 350:0.95, 500:0.84}
    # image_mJy *= conv_factor[band]

    image[0].data = image_mJy

    # update header
    image[0].header['BUNIT'] = 'mJy/pixel'
    image[0].header['history'] = 'units converted to  mJy/pixel'

    return image


def muse_image_to_W_per_m2(image):
    """Converts the units of a muse image (datacube integrated over wavelength)
    lines to W/m^2
    """
    # muse units 10^-20 erg/s/cm^2/A
    # we already integrated in lambda when we did a filter for a given lines
    data = image[0].data
    data = data * (u.erg/u.s/u.cm**2).to(u.W/u.m**2) * 1e-20
    image[0].data = data
    image[0].header['BUNIT'] = 'W/m**2'
    image[0].header['history'] = 'units converted to  W/m**2'
    return image

# other functions


def zero_to_nan(image, value=0):
    """Takes an image whose NaN values had been converted to a number, for
    instance after running a pyraf task, and converts them back to NaN.

    Input:
        name: Name of the image to be fixed.
        value: Value to wich the NaN had been converted to. Default 0.

    Output:
        Image: Image with values equal to value converted to NaN
    """
    data = image[0].data
    w =  np.where(data == value)
    data[w] = np.NaN
    image[0].data = data
    return image

# Functions to run main


def load_fits(name):
    """ Open a fits file image
    Inputs:
        name: name of the .fits file (str).
    Output:
        image:
    """
    while True:
        try:
            file = fits.open(name)
            image = file.copy()
            return image, name
        except FileNotFoundError:
            print(f"File {name} not found")
            name = input('Please enter a different file name: ')


def run_task(image, number, functions):
    """run the function selected by the user. Also verify the user input
    """
    # list of functions with 2 args (image and band)
    f_2_args = [galex_to_mjy_per_pix,
                galex_to_mjy_per_pix_2,
                f_2mass_to_mj_per_pix,
                spitzer_to_mJy_per_pix,
                herschel_to_mJy_per_pix]

    while True:
        try:
            number = int(number)
            break
        except ValueError:
            print('Input not a number')
            print('Choose the number of task to be perform: ' +
                  ', '.join(functions.keys()))
            number = input(': ')

    keys = list(functions.keys())
    keys.sort()
    f = functions[keys[number - 1]]

    if f in f_2_args:
        band = input('Please input the name of the band: ')
        new_image = f(image, band)
    else:
        new_image = f(image)
    return new_image


if __name__ == "__main__":

    # list of functions to choose
    functions = {'1: GALEX to mJy/pix': galex_to_mjy_per_pix_2,
                 '2: 2MASS to mJy/pix': f_2mass_to_mj_per_pix,
                 '3: Spitzer to mJy/pix': spitzer_to_mJy_per_pix,
                 '4: Herschel to mJy/pix': herschel_to_mJy_per_pix,
                 '5: MUSE integrated image to W/m**2': muse_image_to_W_per_m2,
                 '6: CTIO 0.9m CFCCD to mJy/pix': cfccd_to_mjy_per_pix,
                 '7: Zero to nan': zero_to_nan}

    image_name = input('Please enter image name: ')
    image, image_name = load_fits(image_name)

    print('Choose the number of task to be perform: ' +
          ', '.join(functions.keys()))
    task = input(': ')
    new_image = run_task(image, task, functions)

    # to name new image
    sufix = '_converted.fits'
    if task == 6:
        sufix = '_nan.fits'
    name = image_name[:image_name.index('.fits')] + sufix

    # save to file, carefull overwrite set to True
    new_image.writeto(name, overwrite=True)
    image.close()
    print('image saved as: ' + name)
