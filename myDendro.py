from astrodendro import Dendrogram, ppv_catalog
from astrodendro.analysis import PPVStatistic
from spectral_cube import SpectralCube

from astropy.io import fits
from astropy.wcs import WCS


import astropy.units as u
#first step before 

#because this has to be done on server, make this simple

class DendroClass:


	dendroData=None




	def __init__(self,CO12FITS ,dendroFITS,regionName): 
		
		self.CO12FITS=CO12FITS
		
		hdu=fits.open( self.CO12FITS)[0]
		
		self.CO12Data= hdu.data
		
		self.CO12Head= hdu.header
		
		self.dendroFITS=dendroFITS

		self.regionName=regionName


		self.dendroCat=self.regionName+"_DendroCat.fit"


	
	def readDendro(self):

		self.dendroData= Dendrogram.load_from(self.dendroFITS )

	def WriteCatalog(self):
		"""
		
				
		"""
		
		if self.dendroData == None:
			
			self.readDendro()
		
		
		#define metadata
		metadata = {}
		metadata['data_unit'] = u.Kelvin   #u.Jy / u.beam 
		metadata['spatial_scale'] =  30 * u.arcsec
		metadata['beam_major'] =  50 * u.arcsec # FWHM
		metadata['beam_minor'] =  50 * u.arcsec # FWHM
		metadata['wcs'] = WCS( self.CO12Head)  # 22.9 * u.arcsec # FWHM
		
		c= 299792458.
		f=115271202000.0
		wavelength=c/f*u.meter
 
		metadata['wavelength'] = wavelength  # 22.9 * u.arcsec # FWHM

		
		
		cat = ppv_catalog(self.dendroData, metadata)
		
		#write catalog
		try:
			os.remove(self.dendroCat )
		cat.write(self.dendroCat ,format='fits')
		
doDendro= DendroClass( "G130150merge12.fits", "G130150Dendro.fits" "G130150" ) 