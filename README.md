# utils
Python respository for scripts with common fuctions used to manipulate astronomical images. New functions and script will be added as I find myself running the same tasks over and over. 

Scripts:
  Units_conv.py: 
    This script contains functions to convert the units of astronomical images. So far the conversions include:
      - GALEX to mJy/pix (converting to f_lambda first and then to f_nu)
      - Galex to mJy/pix (using zero point values)
      - 2MASS to mJy/pix
      - Spitzer to mJy/pix
      - Herschel to mJy/pix
      - MUSE image (Datacube integrated in wavelength) to W/m^2
    More conversions will be added as needed
    Also includes a function to convert values to NaN, to use in case a routine converts them to zero 
