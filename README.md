# MUSEpy: a python-based analysis tool for MUSE data


This is an example of MUSEpy commands.
Starting from a fully reduced MUSE datacube (DATA extension), obtain [SII]6731, H_alpha and H_beta integrated intensity and velocity maps, compute their line ratio, and make nice figures from these.


import musepy as mp

Step 1: make integrated intensity maps:

mp.moments('cube.fits', [4861.33, 6562.8, 6730.85], ['Hb', 'Ha', 'SII'], 0, save_file=True)


Step 2: make velocity maps, using vacuum wavelengths for the conversion:

mp.moments('cube.fits', [4862.69, 6564.61, 6732.71], ['Hb', 'Ha', 'SII'], 1, save_file=True)


Step 3: make line ratio map:

mp.two_line_ratio('Ha_moment0.fits', 'Hb_moment0.fits', save_fits=True, filename='Ha_Hb.fits')
mp.two_line_ratio('SII_moment0.fits', 'Ha_moment0.fits', save_fits=True, filename='SII_Ha.fits')


Step 4: make pretty pictures:

mp.figure('Ha_moment0.fits', save_file=True, output_file='Ha_moment0.png', minmax=False)
mp.figure('Ha_moment1.fits', save_file=True, output_file='Ha_moment1.png', minmax=False)
mp.figure('SII_moment0.fits', save_file=True, output_file='SII_moment0.png', minmax=False)
mp.figure('SII_moment1.fits', save_file=True, output_file='SII_moment1.png', minmax=True, vmin=-50, vmax=60)
mp.figure('Ha_Hb.fits', save_file=True, output_file='Ha_Hb.png', minmax=False)  
mp.figure('SII_Ha.fits', save_file=True, output_file='SII_Ha.png', minmax=False)



Creation of a subcube based on a DS9 region file:

mp.subcube('cube.fits', 'my_region.reg', save_file=True, filename='subcube.fits')


Compute MUSE telescope time (including overheads), with the number of pointings (NPT), the exposure time per pointing in seconds(DIT) and the number of exposures per pointing (NDIT):

mp.overheads(NPT=9, DIT=90, NDIT=3)

