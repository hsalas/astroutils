import warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from photutils import CircularAperture, CircularAnnulus
# from photutils import SkyCircularAperture, SkyCircularAnnulus
from photutils import centroid_sources, aperture_photometry
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyUserWarning


def ap_an_phot(image, sources, source_ap, sky_ap, delta_ann=3.0, coords='xy',
            centroid=False, show_plot=False, **kargs):
    """Performs circular aperture fotometry on a image, given a table with the
    coordinates (x,y) or (ra,dec) of the sources. Background is estimated from a
    annular aperture and substracted

    Input:
        image:      (str) name of the .fits file with the image
        sources:    table with the Positions of the sources in the image.
        source_ap:  Radius of the circular aperture in pixels
        sky_ap:     inner radius of the annulus to calculate the sky.
        delta_ann:  Width of the annulus to calculate the sky
        centroid:   (Boolean) if True the source position centroid are
                    recalculated.
        sky_coord   (str) type of coordinates to use, either xy or radec.
                    If xy sources table must have 'x' and 'y' cols.
                    If ra dec sources table must to have 'ra' and 'dec' cols.
        show_plot:  (Boolean) If True a plot of the aperture is displayed.
    Output:
        phot_table: Table with the resulting photometry
    """

    data = fits.getdata(image)
    x, y = sources.colnames

    if coords == 'radec':
        sources = SkyCoord(sources[x], sources[y], **kargs)
        # convert to pixel coordinates
        header = fits.getheader(image)
        wcs = WCS(header)
        sources = Table(sources.to_pixel(wcs), names=[x, y])
    elif coords == 'xy':
        warnings.warn('Image pixels start at 1 while python index stats at 0',
                      AstropyUserWarning)
        # Check if the following correction is necessary
        # sources[x] -= 1
        # sources[y] -= 1

    if centroid:
        sources = centroid_sources(data, sources[x], sources[y], box_size=30)
        sources = Table(sources, names=[x, y])

    # create apertures
    sources = [(a, b) for a, b in zip(sources[x], sources[y])]
    source_aperture = CircularAperture(sources, r=source_ap)

    sky_aperture = CircularAnnulus(sources, r_in=sky_ap, r_out=sky_ap +
                                   delta_ann)

    if show_plot:
        index = [str(i+1) for i in range(len(sources))]
        norm = simple_norm(data, 'sqrt', percent=99)
        plt.figure()


        source_aperture.plot(color='white', lw=2)
        sky_aperture.plot(color='red', lw=2)
        for a, b in sources:
            plt.text(a, b, index, color="purple", fontsize=12)
        plt.show()

    # get background using sigma clipped median
    annulus_masks = sky_aperture.to_mask(method='center')
    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
    bkg_median = np.array(bkg_median)

    phot_table = aperture_photometry(data, source_aperture)
    phot_table['annulus_median'] = bkg_median
    phot_table['aper_bkg'] = bkg_median * source_aperture.area()
    phot_table['aper_sum_bkgsub'] = phot_table['aperture_sum'] - \
        phot_table['aper_bkg']

    return phot_table


def ap_phot(image, sources, source_ap, coords='xy', centroid=False,
            show_plot=False, **kargs):
    """Performs circular aperture fotometry on a background substracted image, given a table with the
    coordinates (x,y) or (ra,dec) of the sources.

    Input:
        image:      (str) name of the .fits file with the image
        sources:    table with the Positions of the sources in the image.
        delta_ann:  Width of the annulus to calculate the sky
        centroid:   (Boolean) if True the source position centroid are
                    recalculated.
        sky_coord   (str) type of coordinates to use, either xy or radec.
                    If xy sources table must have 'x' and 'y' cols.
                    If ra dec sources table must to have 'ra' and 'dec' cols.
        show_plot:  (Boolean) If True a plot of the aperture is displayed.
    Output:
        phot_table: Table with the resulting photometry
    """

    data = fits.getdata(image)
    x, y = sources.colnames

    if coords == 'radec':
        sources = SkyCoord(sources[x], sources[y], **kargs)
        # convert to pixel coordinates
        header = fits.getheader(image)
        wcs = WCS(header)
        sources = Table(sources.to_pixel(wcs), names=[x, y])
    elif coords == 'xy':
        warnings.warn('Image pixels start at 1 while python index stats at 0',
                      AstropyUserWarning)
        # Check if the following correction is necessary
        # sources[x] -= 1
        # sources[y] -= 1

    if centroid:
        sources = centroid_sources(data, sources[x], sources[y], box_size=30)
        sources = Table(sources, names=[x, y])

    # create apertures
    sources = [(a, b) for a, b in zip(sources[x], sources[y])]
    source_aperture = CircularAperture(sources, r=source_ap)


    if show_plot:
        index = [str(i+1) for i in range(len(sources))]
        norm = simple_norm(data, 'sqrt', percent=99)
        plt.figure()


        source_aperture.plot(color='white', lw=2)
        for a, b in sources:
            plt.text(a, b, index, color="purple", fontsize=12)
        plt.show()

    phot_table = aperture_photometry(data, source_aperture)

    return phot_table


def bkg_subs(image, box, snr=2, npixels=5, dilate_size=15):
    """Substrackt the background of an image. Sources are idenfified first and
    masked, then the background is estimated as the (sigma-clipped) median.
    """
    # get data from file
    data = fits.getdata(image)
    # create source mask
    mask = make_source_mask(data, snr=snr, npixels=npixels,
                            dilate_size=dilate_size)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, box, filter_size=1, bkg_estimator=bkg_estimator,
                       mask=mask)
    return data - bkg.Background
