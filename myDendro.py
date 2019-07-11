from astrodendro import Dendrogram, ppv_catalog
from astrodendro.analysis import PPVStatistic
from spectral_cube import SpectralCube

from astropy.io import fits
from astropy.wcs import WCS
import os
from astropy.table import Table,vstack
import astropy.units as u
#first step before 

import numpy as np



#because this has to be done on server, make this simple

# oridyce dendrofits

# produce tree , catalog, and save mask files

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

		self.catWithLevelTB=regionName+ self.catWithLevelTB
		self.dendroCat=self.regionName+"_DendroCat.fit"

		self.maskPath=self.regionName+"Mask/"

		

		os.system( "mkdir " +self.maskPath )
		#self.maskPath="./dendroMask/"

	def readDendro(self):

		self.dendroData= Dendrogram.load_from(self.dendroFITS )
		print ""
		print "Dendrogram read!!"
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
		
		if self.regionName=="Mopra":
			metadata['beam_major'] =  33 * u.arcsec # FWHM
			metadata['beam_minor'] =  33 * u.arcsec # FWHM
		
		
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
		
		print "Doing ",cloudFeature
		
		saveFITS=self.maskPath+cloudFeature+"int.fits"
		saveFITSMask=self.maskPath+cloudFeature+"mask.fits"

		

		cloudMask=cloud.get_mask()
		cloudMaskInt=cloudMask.astype('short') #convert mask to int, essentially 0, 1
		intMask= np.sum(cloudMaskInt,axis=0) 
		intMask[intMask>1]=1


		if np.sum(intMask )< 3600 and  cloudLevel <0:
			
			print "Independent small cloud, ignore"
			
			return 

		
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
		
		
		print "saving masks...."
		
		mask_hdu.writeto(saveFITSMask)
		data_hdu.writeto(saveFITS)
	

		

	def produceMaskFITS(self):
		"""
		
		produce mask, in 
		
		"""
		#read cat
		

		if not os.path.isfile(self.dendroCat):
			self.WriteCatalog()
 			
			

		catTB=Table.read(self.dendroCat)
		
		newCatTB=catTB.copy()
		
		newLevelRow= newCatTB.columns[0].copy()
		
		newLevelRow.name="level"
		
		newCatTB.add_column(  newLevelRow )
		
		
		
		fitsHdu=fits.open(self.CO12FITS)[0]
		
 		data=fitsHdu.data
 		head=fitsHdu.header
		
		indexCol= list( newCatTB["_idx"]  )
		
		
		print "Producing mask fits"
		
		if self.dendroData == None:
			
			self.readDendro()
		
		
		for eachC in self.dendroData:
			
			#tempTB=newCatTB[  newCatTB["_idx"]==eachC.idx  ]
			
			rowIndex=indexCol.index(eachC.idx )
			newCatTB[rowIndex]["level"]=eachC.level
			
 
			self.saveMask( eachC,data,head)
		#for eachC in self.dendroData:

 

		try:
			os.remove(self.catWithLevelTB)
		except:
			pass

		newCatTB.write(self.catWithLevelTB )



	def calDisByCloudID(self):
		"""
		"""

	def writeTreeStructure(self):
		
		if self.dendroData == None:
			
			self.readDendro()
		
		f=open("treeStructure.txt",'w')
		
		#for eachC in self.dendroData:
		for eachC in self.dendroData:

			parentID=-1
			
			p=eachC.parent
			
			if p!=None:
			
				parentID=p.idx
			
			fileRow="{} {}".format(eachC.idx,parentID)
			f.write(fileRow+" \n")
			
		f.close()

	def ZZZ(self):
		pass


if 1:
	# this should only be running on the server
	doDendro= DendroClass( "/home/qzyan/WORK/dataDisk/G190/G190MergeCO12_32bit.fits", "G190Dendro32bit.fits","G190" ) 
	
	#doDendro.WriteCatalog()
	
	doDendro.produceMaskFITS()
	



if 0:
	doDendro= DendroClass( "G130150merge12.fits", "G130150Dendro.fits","G130150" ) 
	doDendro.writeTreeStructure()

if 0:
	# this should only be running on the server
	doDendro= DendroClass( "G130150merge12.fits", "G130150Dendro.fits","G130150" ) 
	
	#doDendro.WriteCatalog()
	
	doDendro.produceMaskFITS()

 
 
if 0:
	# this should only be running on the server
	doDendro= DendroClass( "moproCO12.fits", "MopraDendroLarge.fits","Mopra" ) 
	
	#doDendro.WriteCatalog()
	
	doDendro.produceMaskFITS()
	