from astropy.io import fits
from astropy.visualization import simple_norm
import aplpy

def load_fits(name):
    """ Function to open a fits file a return its data and header
    Inputs:
        name: name of the .fits file (str).
    Output:
        data: image data
        hdr:  image header
        name: image name
    """
    while True:
        try:
            # image = fits.open(name)
            data = fits.getdata(name)
            hdr = fits.getheader(name)
            return data, hdr, name
        except FileNotFoundError:
            print(f"File {name} not found")
            name = input('Please enter a different file name: ')


def find_pixel_scale(header):
    """Finds the value of the image pixel scale from the image headers
    Inputs:
        header: Header of the image
    Output:
        pixel_scale: Pixel scale of the image in arcsec/pixel
    """

    pixel_scale = None
    keys = [key for key in header.keys()]
    if ('CD1_1' in keys) and ('CD1_2' in keys):
        pixel_scale = np.sqrt(header['CD1_1']**2 + header['CD1_2']**2)*3600

    elif ('PC1_1' in keys) and ('PC2_2' in keys):
        pixel_scale = np.sqrt(header['PC1_1']**2 + header['PC2_2']**2)*3600

    elif 'PXSCAL_1' in keys:
        pixel_scale = abs(header['PXSCAL_1'])

    elif 'PIXSCALE' in keys:
        pixel_scale = header['PIXSCALE']

    elif 'SECPIX' in keys:
        pixel_scale = header['SECPIX']

    elif 'CDELT1' in keys:
        pixel_scale = abs(header['CDELT1'])*3600

    else:
        print('Unable to get pixel scale from image header')
        while True:
            pixel_scale = input('Plesae input the pixel scale value in \
                                arcsec per pixel')
            try:
                pixel_scale = float(pixel_scale)
                return pixel_scale
            except ValueError:
                pass

    return pixel_scale


def plot_fits(name, scale=None, cmap='default', colorbar='none',
              **kwargs):
    """ plots a intensity map in (ra, dec) of a fits image.

    Inputs:
        name:   Name of the fitsfile (str)
        scale:  Scalebar marker to add to the plot. Tuple (d,t) with d the
                distance (float) and t (str) the text to display
        cmap:   Colormap to be used (str). Check list of colormps in
                available in matplotlib https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
        colorbar:   Position of the colorbar (str), options: 'top', 'bottom',
                    'left', 'right' or 'none'(default)
        **kwargs:   kargs from aplpy and matplotlib should be recognized.

    Returns:
        image plot
    """
    with fits.open(name) as f:
        data = f[0].data
    maxval = max(np.abs(np.nanmin(data)), np.abs(np.nanmax(data)))
    fig = aplpy.FITSFigure(name, **kwargs)
    fig.show_colorscale(vmin=-maxval, vmax=maxval, cmap=cmap)
    if scale != None:
        fig.add_scalebar(scale[0])
        fig.scalebar.set_label(scale[1])
        fig.scalebar.set_color('white')
        fig.set_nan_color('grey')
    if colorbar =='top':
        fig.add_colorbar()
        fig.colorbar.set_location('top')
    elif colorbar =='bottom':
        fig.add_colorbar()
        fig.colorbar.set_location('bottom')
    elif colorbar =='left':
        fig.add_colorbar()
        fig.colorbar.set_location('left')
    elif colorbar =='right':
        fig.add_colorbar()
        fig.colorbar.set_location('right')
    # fig.set_theme('publication')
    fig.frame.ax.figure.tight_layout()
    return fig


def get_idx(hdr, l, side='left'):
    """ Returns the index of a corresponding to wavelength 'l' in a spectra.

    Inputs:
        hdr:    Header of the .fits file with the spectra.
        l:      the desired wavelength.
        side :  {'left', 'right'}, optional
                If 'left', the index of the first suitable location found is
                given. If 'right', return the last such index. If there is no
                suitableindex, return either 0 or N
    Output:
        idx:    indes corresponding to wavelength l
    """
    #read from header

    keys = [key for key in hdr.keys()]

    ref_p = hdr['CRPIX3']
    ref_l = hdr['CRVAL3']
    n_pix = hdr['NAXIS3']

    if 'CD3_3' in keys:
        dl = hdr['CD3_3']
    elif 'CDELT3' in keys:
        dl = hdr['CDELT3']
    else:
        raise ValueError('keywords not found in header')

    if 'PC3_3' in keys:
        conv = hdr['PC3_3']
        dl *= conv

    ini = ref_l - (ref_p - 1) * dl
    fin = ref_l + (n_pix - ref_p) * dl
    aux = np.arange(ini, fin, dl)
    idx = np.searchsorted(aux, l, side=side)

    return idx
