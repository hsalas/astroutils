from astropy.io import fits

def load_fits(name):
    """ Function to open a fits file a return its data and header
    Inputs:
        name: name of the .fits file (str).
    Output:
        data: image data
        hdr:  image header
    """
    while True:
        try:
            # image = fits.open(name)
            data = fits.getdata(name)
            hdr = fits.getheader(name)
            return data, hdr
        except FileNotFoundError:
            print(f"File {name} not found")
            name = input('Please enter a different file name: ')
