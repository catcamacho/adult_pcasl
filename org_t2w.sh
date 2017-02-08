#! /bin/csh

set aslpath = /share/iang/active/BABIES/BABIES_perfusion/subjDir/rawasl
set cbfpath = /share/iang/active/BABIES/BABIES_perfusion/subjDir/postproc
set ibeatpath = /share/iang/active/BABIES/BABIES_ibeat/subjDir
set log = /share/iang/active/BABIES/BABIES_perfusion/subjDir/orgLog.txt

foreach sub(002x 012 020 021 027 031 032 033x 035 036 045)
set folder = ${sub}-BABIES-T1
set ibeat = T1${sub}
if (-e $ibeatpath/$ibeat/${ibeat}-5/skullstripped_anat.nii) then
	echo "----------- copying anat for subject " $sub >> $log

	cp $ibeatpath/$ibeat/${ibeat}-5/skullstripped_anat.nii $cbfpath/$folder/t2anat.nii
	cp $ibeatpath/$ibeat/${ibeat}-5/skullstripped_anat.nii $aslpath/$folder/t2anat.nii
	
	echo "----------- Subject " $sub " is copied!" >> $log
else
	echo "WARNING: No anat found for subject " $sub >> $log
endif

end
	
	