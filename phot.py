from astropy.io import fits
from astropy.table import Table
from photutils import centroid_sources, CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm

def ap_phot(image, sources_pos, source_ap, sky_ap, delta_ann=3.0,
            centroid=False):
    """Performs circular aperture fotometry on a image, given the x,y position of the sources.

    Input:
        image: Image data.
        sources_pos: Positions of the sources in the image.
        source_ap:  Radius of the circular aperture.
        sky_ap:     inner radius of the annulus to calculate the sky.
        delta_ann:  Width of the annulus to calculate the sky
        centroid:   (Bolean) if True the source position centroid are
                    recalculated.
    Output:
        phot_table: Table with the resulting photometry
    """
    # read list of sources
    sources = Table.read(sources_pos, format='ascii', names=['x','y'])

    # python index starts at zero
    sources['x'] -= 1
    sources['y'] -=1

    if centroid:
        sources = centroid_sources(image, sources['x'], sources['y'])

    # create apertures
    source_aperture = CircularAperture((sources['x'],sources['y']), r =
                                       source_ap)

    sky_aperture = CircularAnnulus((sources['x'],sources['y']), r_in = sky_ap
                                   , r_out = sky_ap + delta_ann)

    # get background using sigma clipped median
    annulus_masks = sky_aperture.to_mask(method='center')
    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(image)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
    bkg_median = np.array(bkg_median)


    phot_table = aperture_photometry(image, source_aperture)
    phot_table['annulus_median'] = bkg_median
    phot_table['aper_bkg'] = bkg_median * source_aperture.area()
    phot_table['aper_sum_bkgsub'] = phot_table['aperture_sum'] - \
                                    phot_table['aper_bkg']

    # norm = simple_norm(image, 'sqrt', percent=99)
    # plt.imshow(image, norm=norm)
    # source_aperture.plot(color='white', lw=2)
    # sky_aperture.plot(color='red', lw=2)

    return phot_table
