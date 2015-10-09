"""
MUSEpy is a python based package for the analysis of MUSE data.

Please report requests/bugs/issues to amcleod@eso.org.

Please aknowledge MUSEpy if you use it for your publications!
"""

class MUSEpy:
    def __init__(self):
        self
        

from astropy.io import fits
from spectral_cube import SpectralCube
from astropy import units as u
import aplpy


def figure(image, colorbar=True, colorbar_text='', save_file=False, output_file='out.png', minmax=False, vmin='', vmax=''):
    """

    
    image: str
        The image in fits format to open.

    colorbar: bool, optional
        If True, a colorbar will be added. Default is True.
    
    colorbar_text: str, optional
        The preferred label of the colorbar, only works if colorbar=True. Default is an empty string.

    save_file: bool, optional    
        If True, the figure will be saved. Default is False.

    output_file: str, optional
        The name of the output file, specifying the filetype. Default is 'out.png'.

    minmax: bool, optional
        Manually set minimum and maximum values. If False, minimum and maximum will be set automatically.

    vmin, vmax: float, optional
        Minimum and maximum values

    example:
        image = aplpy_figure(image='image.fits', colorbar_text='text for my colorbar', save_file=True, output_file='output.png')


    
    """

    
    im=aplpy.FITSFigure(image)
    im.show_colorscale()
    if minmax==True:
        im.show_colorscale(vmin=vmin, vmax=vmax)
    if colorbar==True:
        im.show_colorbar()
        im.colorbar.set_axis_label_text(colorbar_text)
    if save_file==True:
        im.save(output_file)
        
    return im


def moments(cube_fits, line_values, line_names, moment, save_file=False):
    """
    cubes: str
        The datacube in fits format to open

    line_values: list of floats
        The wavelengths of the lines. !!! In general: if moment=0 the required wavelenth should be air, if moment=1 it should be vacuum !!!

    line_names: list of str
        The identifier of the lines

    moment: 0 or 1

    save_file: bool, optional
        Set to True if the result is to be saved as a fits file. Default is False.

    example:
        moment = moments('cube.fits', [4861.33, 6562.8], ['Hb', 'Ha'], moment=0)
    """


    print line_values, line_names
    cube=SpectralCube.read(cube_fits)

    for line,stri in zip(line_values,line_names):
        if moment==0:
            mom = cube.spectral_slab((line-3)*u.AA, (line+3)*u.AA).sum(axis=0)
            if save_file==True:
                mom.hdu.writeto(str(stri)+'_moment0.fits',clobber=True)
        if moment==1:    
            mom = cube.with_spectral_unit(u.km/u.s, rest_value=line*u.AA,velocity_convention='optical').spectral_slab(-300*u.km/u.s,300*u.km/u.s).moment1()
            if save_file==True:
                mom.hdu.writeto(str(stri)+'_moment1.fits',clobber=True)

    return mom


def two_line_ratio(line1, line2, save_fits=False, filename=''):
    """
    line1: str
        Filename of nominator line

    line2: str
        Filename of denominator line

    filename: str
        The name of the fits file to be saved 

    Warning: works only if line1 and line2 have the same shape!

    Example:
        mp.two_line_ratio('Ha_moment0.fits', 'Hb_moment0.fits', save_fits=True, filename='Ha_Hb.fits')
    
    """

    
    l1, l2 = fits.getdata(line1), fits.getdata(line2)

    hd = fits.getheader(line1)
    ratio = l1/l2

    if save_fits==True:
        ff = fits.PrimaryHDU(data=ratio, header=hd)
        ff.writeto(filename, clobber=True)

    return ratio

def three_line_ratio(line1, line2, line3, header, savefits, filename, mean):
    """
    Ratio of the type (L1+L2)/L3 if mean=False
    or of the type L12/L3 if mean=True, where L12 = (L1+L2)/2
    
    line1: str or moment object
        Filename of nominator line 1

    line2: str or moment object
        Filename of denominator line 2

    line2: str or moment object
        Filename of denominator line 

    header: .txt file
        The header to be used when saving the ratio to a fits file.
        If not specified, the header of line 1 will be used.    

    filename: str
        The name of the fits file to be saved 

    Warning: works only if line1, line2 and line3 have the same shape!
    
    """
    l1, l2 = fits.getdata(line1), fits.getdata(line2)
    hd = fits.getheader(line1)
    if mean==True:
        ratio = ((l1+l2)/2.)/l3
    if mean==False:
        ratio = (l1+l2)/l3

    if savefits==True:
        ff = fits.PrimaryHDU(data=ratio, header=hd)
        ff.writeto(filename, clobber=True)

    return ratio




def overheads(NPT, DIT, NDIT):
    """

    Compute total telescope time (including instrument overheads)

    NPT = number of poinitings
    DIT = exposure time
    NDIT = number of exposures per pointing

    """
    ov = 360. + 120. + NPT*NDIT*(DIT + 80. + 15.)
    print 'Telescope time in h = ', ov/3600.
