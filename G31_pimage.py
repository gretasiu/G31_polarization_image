import os
import numpy as np

from astropy.io.fits import getdata
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import constants as con
from astropy.coordinates import SkyCoord

import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib

import aplpy

matplotlib.use('PDF')

intensitymap = 'G31p4_Qband_D.rob2.I.image.tt0.dropdeg.fits'
intensity_scale = 1.0
if_success = False
try:

    # importing FITS image to a HDU
    Ihdu = fits.open(intensitymap)

    # editing the FITS image by multiplying a scaling factor
    Ihdu[0].data = Ihdu[0].data * intensity_scale

    if_success = True

except:
    print('Unable to read the intensity FITS image. Please double check the image file.')
    print(intensitymap)

if (if_success == True):
    # Reading FITS header
    try:
        naxis1 = Ihdu[0].header['naxis1']
        naxis2 = Ihdu[0].header['naxis2']
        crval1 = Ihdu[0].header['crval1']
        crpix1 = Ihdu[0].header['crpix1']
        cdelt1 = Ihdu[0].header['cdelt1']
        crval2 = Ihdu[0].header['crval2']
        crpix2 = Ihdu[0].header['crpix2']
        cdelt2 = Ihdu[0].header['cdelt2']
        hduwcs = wcs.WCS(Ihdu[0].header)
    except:
        print('Warning. No coordinate headers')

    try:
        bmaj = Ihdu[0].header['bmaj']
        bmin = Ihdu[0].header['bmin']
        bpa = Ihdu[0].header['bpa']
    except:
        print('Warnning. No header for synthesized beam size')

# pixel values are in Ihdu[0].data[x][y]
# value = Ihdu[0].data[240][240]
# print(value)
# Ihdu .info()


def plot_intensity(figsize=[9.0, 9.0],
                   plot_intensity=True,
                   intensity_scale=1.0,
                   vmax=-999.0, vmin=0.0,
                   cmap='viridis',
                   plot_ticks=True,
                   tick_font=15,
                   tick_ypad=0,
                   plot_colorbar=True,
                   colorbar_location='top',
                   colorbar_width=0.2, colorbar_pad=0.1, colorbar_font=15,
                   colorbar_label='Colorbar', colorbar_labelpad=12, colorbar_labelfont=15,
                   ra_center=1, dec_center=1, width=1.0, height=1.0,
                   plot_scalebar=False,
                   distance=140.0,
                   scalebar_size=100.0, scalebar_text='100 au',
                   scalebar_color=(0, 0, 0),
                   scalebar_font=25.0, scalebar_linewidth=3.0,
                   plot_beam=True, beam_color='black'
                   ):

    if_plot = False
    try:
        fig = aplpy.FITSFigure(Ihdu, figsize=(figsize[0], figsize[1]))
        if_plot = True

    except:
        print('Unable to load or plot hdu. Please double check the image file.')

    if (if_plot == True):
        if (vmax == -999.0):
            vmax = np.nanmax(Ihdu[0].data)*1.2
            vmin = 0.0

        fig.show_colorscale(
            vmax=vmax, vmin=vmin,
            interpolation='bicubic',
            cmap=cmap
        )

    # ticks
    fig.tick_labels.set_font(size=tick_font)
    fig.axis_labels.set_font(size=tick_font)
    fig.axis_labels.set_ypad(tick_ypad)

    if (plot_ticks != True):
        fig.axis_labels.hide_x()
        fig.axis_labels.hide_y()
        fig.tick_labels.hide_x()
        fig.tick_labels.hide_y()

    # recentering
    if ((ra_center != 0) or (dec_center != 0)):
        print('recentering')
        #fig.recenter(ra_center, dec_center, width=width, height=height)
        fig.recenter(240, 240, width=0.3, height=0.2)  # degrees

    # plot color bar
    if (plot_colorbar == True):
        fig.add_colorbar()
        fig.colorbar.show()
        fig.colorbar.set_location(colorbar_location)
        fig.colorbar.set_width(colorbar_width)
        fig.colorbar.set_pad(colorbar_pad)
        fig.colorbar.set_font(size=colorbar_font, weight='medium',
                              stretch='normal', family='sans-serif',
                              style='normal', variant='normal')

        fig.colorbar.set_axis_label_text(colorbar_label)
        fig.colorbar.set_axis_label_font(size=colorbar_labelfont)
        fig.colorbar.set_axis_label_pad(colorbar_labelpad)

    # plot synthesized beam
    if (plot_beam == True):
        fig.add_beam()
        fig.beam.set_color(beam_color)

    # plot scalebar
    if (plot_scalebar == True):

        fig.add_scalebar(scalebar_size * (1.0/distance) / 3600.0)
        fig.scalebar.set_label(scalebar_text)
        fig.scalebar.set_color(scalebar_color)
        fig.scalebar.set_font(size=scalebar_font)
        fig.scalebar.set_linewidth(scalebar_linewidth)

    outfigname = 'G31_I.pdf'
    fig.save(outfigname, transparent=True)


plot_intensity(colorbar_label='(Jy/beam)')
