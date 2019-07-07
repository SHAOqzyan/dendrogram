from astrodendro import Dendrogram, ppv_catalog
from astrodendro.analysis import PPVStatistic
from spectral_cube import SpectralCube

from astropy.io import fits
from astropy.wcs import WCS
import os

import astropy.units as u
#first step before 

#because this has to be done on server, make this simple

class DendroClass:


	dendroData=None

	CO12ResV= 0.158737644553 #km/s


	catWithLevelTB="cloudCatWithLevel.fit"
	
	def __init__(self, CO12FITS ,dendroFITS,regionName): 
		
		self.CO12FITS=CO12FITS
		
		hdu=fits.open( self.CO12FITS)[0]
		
		self.CO12Data= hdu.data
		
		self.CO12Head= hdu.header
		
		self.dendroFITS=dendroFITS

		self.regionName=regionName


		self.dendroCat=self.regionName+"_DendroCat.fit"

		self.maskPath="dendroMask"

		os.system( "mkdir " +self.maskPath )
		self.maskPath="./dendroMask/"

	def readDendro(self):

		self.dendroData= Dendrogram.load_from(self.dendroFITS )

	def WriteCatalog(self):
		"""
		
				
		"""
		
		print "Writing catalg"
		
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
		except:
			pass
		cat.write(self.dendroCat ,format='fits')

	def headerTo2D(self,head):
		"""
		"""
		newHead=head.copy()
		del newHead["CRPIX3"]
		del newHead["CDELT3"]
		del newHead["CRVAL3"]
		del newHead["CTYPE3"] 
		
		newHead["NAXIS3"]=2
		newHead["WCSAXES"]=2
 
		return newHead

	def saveMask(self,cloud, FITSData,fitsHeader ):
		
		"""
		save Mask fits 
		"""
		cloudID=cloud.idx
		cloudLevel=cloud.level
		
		#saveName=self.dataPath+"dendroCloud{}_level{}.png".format(cloudID,cloudLevel)
		
		cloudFeature="Cloud{}_level{}".format(cloudID,cloudLevel)
		
		
		saveFITS=self.maskPath+cloudFeature+"int.fits"
		saveFITSMask=self.maskPath+cloudFeature+"mask.fits"



		cloudMask=cloud.get_mask()
		cloudMaskInt=cloudMask.astype('short') #convert mask to int, essentially 0, 1
		intMask= np.sum(cloudMaskInt,axis=0) 
		intMask[intMask>1]=1

		
		sumXY= np.sum(cloudMaskInt,axis=(1,2)) 
		
		nonZeros= np.nonzero(sumXY)[0]
		
		first= 	nonZeros[0]
		second= nonZeros[-1] 
		
		
		cloudCut=FITSData[first:second+1,:,:]
		#print first, second
 
		intData=np.sum(cloudCut,axis=0)*self.CO12ResV #convert to K km/s 
		
		intHead=self.headerTo2D( fitsHeader )
		
	 
		mask_hdu = fits.PrimaryHDU(intMask, intHead)

		data_hdu = fits.PrimaryHDU(intData, intHead)
		
		
		mask_hdu.writeto(saveFITSMask)
		data_hdu.writeto(saveFITS)
	

		

	def produceMaskFITS(self):
		"""
		
		produce mask, in 
		
		"""
		#read cat
		
		
		catTB=Table.read(self.dendroCat)
		
		newCatTB=catTB.copy()
		
		newLevelRow= newCatTB.columns[0].copy()
		
		newLevelRow.name="level"
		
		newCatTB.add_column(  newLevelRow )
		
		
		
		fitsHdu=fits.open(searchFITS)[0]
		
 		data=fitsHdu.data
 		head=fitsHdu.header
		
		indexCol= list( newCatTB["_idx"]  )
		
		for eachC in self.dendroData:
			
			#tempTB=newCatTB[  newCatTB["_idx"]==eachC.idx  ]
			
			rowIndex=indexCol.index(eachC.idx )
			newCatTB[rowIndex]["level"]=eachC.level
			
			statsC =PPVStatistic( eachC )
			
			if eachC.level>0 and statsC.area_exact.value<3600:
				continue # dot con 
 	
			#else write the mask fits
			
			#cloudID="Cloud{}".format(eachC.idx)
			
			self.saveMask( eachC,data,head)
		#for eachC in self.dendroData:

 

		try:
			os.remove(self.catWithLevelTB)
		except:
			pass

		newCatTB.write(self.catWithLevelTB )

doDendro= DendroClass( "G130150merge12.fits", "G130150Dendro.fits","G130150" ) 

#doDendro.WriteCatalog()

doDendro.produceMaskFITS()