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
	
	
	G130150Mask="./G130150Mask/"
	
	def __init__(self, CO12FITS ,dendroFITS,regionName): 
		
		self.CO12FITS=CO12FITS
		

		
		self.dendroFITS=dendroFITS
		
		#read fits imediately
		
		

		

		hdu=fits.open( self.CO12FITS)[0]
		
		self.CO12Data= hdu.data
		
		self.CO12Head= hdu.header
		
		self.regionName=regionName

		self.catWithLevelTB=regionName+ self.catWithLevelTB
		self.dendroCat=self.regionName+"_DendroCat.fit"

		self.maskPath=self.regionName+"Mask/"
		os.system( "mkdir " +self.maskPath )
		
		self.pvPath= "./pvPath/"+self.regionName+"PVFITS/"
		os.system( "mkdir " +self.pvPath )


		print "Loading dendrodata..."
		self.dendroData= Dendrogram.load_from(self.dendroFITS )

		#self.maskPath="./dendroMask/"

	def readDendro(self):
		
		print "Reading dendro saved fits"
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
		# the velocity_scale should be mentioned, other wise the v_rms would be in pixels
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

		try:
			os.remove( saveFITSMask )
		except:
			pass

		try:
			os.remove( saveFITS )
		except:
			pass

		mask_hdu.writeto(saveFITSMask)
		data_hdu.writeto(saveFITS)
	
	
	
	
	
	
	
	def producePVFITS(self,CO12FITS,pvHeaderFITS,dendroTB=None   ):
		
		"""
		
		produce PV FITS
		
		the pvHeaderFITS is only used for provide header fits, which could be produced with myPYTHON (creatPPVHeader)
		
		"""
		
		#origin dat
 
		searchFITS= CO12FITS
		#searchFITS13="./data/mosaic_L.fits"
		

		fitsData = self.CO12Data
		fitsHead =  self.CO12Head

 


		hduPV=fits.open(pvHeaderFITS)[0]

		pvData=  hduPV.data
		pvHeader =   hduPV.header
		

 
		

		#print "11111111111"
		for eachC in self.dendroData:
			#print "222222222222"

			#if eachC.level>0: #only trunks #and statsC.area_exact.value<3600: on
				#continue
			
			#statsC =PPVStatistic( eachC )
 
			
			cloudID="Cloud{}".format(eachC.idx)
			#print dir(statsC)
			

			
			maskData=eachC.get_mask()
			maskData=maskData*1
			cloudData12=fitsData*maskData
			
			PVData2D=np.sum(cloudData12,axis=1) 
			
			PVDataMask2D=np.sum(maskData,axis=1)

			PVDataMask2D[PVDataMask2D>0]=1
			
			pvSaveName=self.pvPath+ "{}_{}_PV.fits".format( self.regionName , cloudID)
			
			pvSaveNameMask=self.pvPath+ "{}_{}_PVMask.fits".format(  self.regionName , cloudID)


			
			print "Saving...",cloudID
			
			fits.writeto(pvSaveName, PVData2D, pvHeader, overwrite=True)
			
			fits.writeto(pvSaveNameMask, PVDataMask2D, pvHeader, overwrite=True)

 
			#self.drawDendroPVbyMaskData(eachC.idx,  fitsData, fitsHead , maskData,pvHeader  )

			 


		return 




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
		
		f=open(self.regionName+"treeStructure.txt",'w')
		
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

if 1:# PV fits

	runi=1
	dendroFITSPath=""
	dendroFITS=  "left{}Dendro.fits".format(runi)
	regionName="left{}".format(runi)
	FITSCO12=  "left{}Cut.fits".format(runi)

	pvHeadFITS= "Gleft{}_LV.fits".format(runi)

	doDendro= DendroClass( FITSCO12,dendroFITS,regionName )

	doDendro.producePVFITS (FITSCO12,pvHeadFITS )


if 0: #intialize a region
	runi=2
	dendroFITSPath=""
	dendroFITS=  "left{}Dendro.fits".format(runi)
	regionName="left{}".format(runi)
	FITSCO12=  "left{}Cut.fits".format(runi)

	doDendro= DendroClass( FITSCO12,dendroFITS,regionName )
	doDendro.WriteCatalog()
	doDendro.produceMaskFITS()
	doDendro.writeTreeStructure()

	#doDendro.producePVFITS(doDendro.CO12FITS, "/home/qzyan/WORK/dataDisk/MWISP/G40/fred12LVHeader.fits" )









if 0: #intialize a region

	dendroFITSPath="/home/qzyan/WORK/projects/doDendro/"
	dendroFITS=dendroFITSPath + "G214Dendro.fits"
	regionName="G214"
	FITSCO12=  "/home/qzyan/WORK/projects/myDendro/G214CO12.fits"

	doDendro= DendroClass( FITSCO12,dendroFITS,regionName )
	#doDendro.producePVFITS(doDendro.CO12FITS, "/home/qzyan/WORK/dataDisk/MWISP/G40/fred12LVHeader.fits" )

	doDendro.produceMaskFITS()

	doDendro.WriteCatalog()
	#doDendro.writeTreeStructure()



if 0: #intialize a region

	dendroFITSPath="/home/qzyan/WORK/projects/doDendro/"
	dendroFITS=dendroFITSPath + "split11G2650dendro.fits"
	regionName="split11New"
	FITSCO12=dendroFITSPath+ "split11G2650.fits"

	doDendro= DendroClass( FITSCO12,dendroFITS,regionName )
	#doDendro.producePVFITS(doDendro.CO12FITS, "/home/qzyan/WORK/dataDisk/MWISP/G40/fred12LVHeader.fits" )

	doDendro.produceMaskFITS()

	doDendro.WriteCatalog()
	doDendro.writeTreeStructure()






if 0: #intialize a region

	dendroFITSPath="/home/qzyan/WORK/projects/doDendro/"
	dendroFITS=dendroFITSPath + "split11G2650dendro.fits"
	regionName="split11New"
	FITSCO12=dendroFITSPath+ "split11G2650.fits"

	doDendro= DendroClass( FITSCO12,dendroFITS,regionName )
	#doDendro.producePVFITS(doDendro.CO12FITS, "/home/qzyan/WORK/dataDisk/MWISP/G40/fred12LVHeader.fits" )

	doDendro.produceMaskFITS()

	doDendro.WriteCatalog()
	doDendro.writeTreeStructure()


if 0: # G40
	doDendro= DendroClass( "fred12.fits", "G40Dendro32bit.fits","G40" )
	doDendro.producePVFITS(doDendro.CO12FITS, "/home/qzyan/WORK/dataDisk/MWISP/G40/fred12LVHeader.fits" )


	#doDendro.produceMaskFITS()

	#doDendro.WriteCatalog()
	#doDendro.writeTreeStructure()


if 0: #producePVFITS
	doDendro= DendroClass( "G130150merge12.fits", "G130150Dendro.fits","G130150" ) 
	doDendro.producePVFITS(doDendro.CO12FITS, "G130150LVHead.fits","cloudCatWithLevelG130150.fit" )


if 0:
	doDendro= DendroClass( "/home/qzyan/WORK/projects/maddalena/data/mosaic_U.fits", "/home/qzyan/WORK/projects/maddalena/G210DendrogramMoreComplete.fits","G210" ) 

	#doDendro.writeTreeStructure()

if 0:
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
	