#! /bin/csh

set origfp = /Volumes/group/iang/active/BABIES/BABIES-T1
set aslpath = /Volumes/group/iang/active/BABIES/BABIES_perfusion/subjDir/rawasl
set cbfpath = /Volumes/group/iang/active/BABIES/BABIES_perfusion/subjDir/postproc
set datapath = functional/asl
set log = /Volumes/group/iang/active/BABIES/BABIES_perfusion/subjDir/orgLog.txt

foreach sub(002x-BABIES-T1 010-BABIES-T1 012-BABIES-T1 020-BABIES-T1 021-BABIES-T1 025-BABIES-T1 027-BABIES-T1 028x-BABIES-T1 031-BABIES-T1 032-BABIES-T1 033x-BABIES-T1 035-BABIES-T1 036-BABIES-T1 040-BABIES-T1 045-BABIES-T1)

if (-e $origfp/$sub/$datapath/cbf) then
	echo "----------- Working on subject " $sub >> $log
	mkdir $aslpath/$sub
	mkdir $cbfpath/$sub
	
	cd $origfp/$sub/$datapath
	
	foreach seq (pcasl pd cbf)
		cd $seq
		to3d -prefix $seq I*
		3dAFNItoNIFTI -prefix $seq ${seq}+orig
		rm ${seq}+orig*
		cd ..
	end
	
	cp $origfp/$sub/$datapath/cbf/cbf.nii $cbfpath/$sub/
	cp $origfp/$sub/$datapath/pcasl/pcasl.nii $aslpath/$sub/
	cp $origfp/$sub/$datapath/pd/pd.nii $aslpath/$sub/
	
	echo "----------- Subject " $sub " is ready!" >> $log
else
	echo "WARNING: No ASL found for subject " $sub >> $log
endif
	
	