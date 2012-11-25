begin
	FileList=new(17,"string")
	FileList( 0)="sfc.nc"
	FileList( 1)="sLAI.nc" 
	FileList( 2)="sNDVI.nc"
	FileList( 3)="sTs.nc"
	FileList( 4)="shc.nc"
	FileList( 5)="sPressure.nc"
	FileList( 6)="sTa_ref.nc"
	FileList( 7)="sqa_ref.nc"
	FileList( 8)="sUref.nc"
	FileList( 9)="sAlbedo.nc"
	FileList(10)="sEmissivity.nc"
	FileList(11)="sSWd.nc"
	FileList(12)="sLWd.nc"
	FileList(13)="shpbl.nc"
	FileList(14)="sSWd24.nc"
	FileList(15)="sTa_avg.nc"
	FileList(16)="sNs.nc"
	
	fids=addfiles(FileList,"r")
	
	fid=addfile("sLSM.nc","r")
	print("getting LandSeaMask")
	LSM=fid->LSM
	ComptMask=where(LSM.ne.1,0,1);
	delete(LSM)
	delete(fid)
	
	print("Loading fc,NDVI, and hc ...")
	fc=fids[0]->fc
	NDVI=fids[2]->NDVI
	hc=fids[4]->hc
	print("Loading LAI ...")
	LAI=fids[1]->LAI
	Albedo=fids[9]->Albedo
	
	print ( "- Applying fc,NDVI,hc, Albedo, and LAI mask ...")
	print(any(fc.eq.-2))
	ComptMask=where(ismissing(fc),0,ComptMask)
	ComptMask=where(ismissing(NDVI),0,ComptMask)
	ComptMask=where(ismissing(hc),0,ComptMask)
	ComptMask=where(ismissing(Albedo),0,ComptMask)
	ComptMask=where(ismissing(LAI),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	print("- Applying Temperatures Mask...")
	tmpVar1=(/ fids[3]->Ts /)
	tmpVar2=(/ fids[6]->Ta_ref /)
	ComptMask=where(ismissing(tmpVar1).or.\
					   ismissing(tmpVar2).or.\
					   tmpVar2.ge.tmpVar1,0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	print("- applying Velocity Mask ...")
	tmpVar1=(/ fids[8]->Uref /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	
	print("- applying Pressure Mask ...")
	tmpVar1=(/ fids[5]->Pressure /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	print("- applying Pressure MSL Mask ...")
	
	print("- applying specific humidity Mask ...")
	tmpVar1=(/ fids[7]->qa_ref /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	
	print("- applying specific Emissivity Mask ...")
	tmpVar1=(/ fids[10]->Emissivity /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	
	print("- applying specific SWd Mask ...")
	tmpVar1=(/ fids[11]->SWd /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	
	print("- applying specific LWd Mask ...")
	tmpVar1=(/ fids[12]->LWd /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	
	print("- applying specific hpbl Mask ...")
	tmpVar1=(/ fids[13]->hpbl /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	
	print("- applying specific SWd24 Mask ...")
	tmpVar1=(/ fids[14]->SWd24 /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	
	print("- applying specific Ta_avg Mask ...")
	tmpVar1=(/ fids[15]->Ta_avg /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	
	print("- applying specific Ns Mask ...")
	tmpVar1=(/ fids[16]->Ns /)
	ComptMask=where(ismissing(tmpVar1),0,ComptMask)
	print("-- Num Cells for Calculation: "+sum(ComptMask)+" => %"+tofloat(sum(ComptMask))/(120.0*120.0)*100.0)	
	
	print(" ");
	
	ComptMask!0="lat"
	ComptMask!1="lon"
	ComptMask@MODISTILE="V04H09"
	ComptMask@Year="2008"
	fid=addfile("sComptMask.nc","c")
	fid->ComptMask=ComptMask
end

