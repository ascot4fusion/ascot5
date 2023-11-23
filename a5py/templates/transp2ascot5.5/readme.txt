The primary script for generating input is **prepare_input.py**


It is assumed that the user has access to TRANSP output files. 

Only three TRANSP files are necessary for this package:

1) The **full CDF file**, of the form "**<run_ID>.CDF**"
2) The "**plasma state file**", of the form "**<run_ID>_ps_ts1_state.cdf**"
3) IF using NUBEAM's markers as input to ASCOT5, the "**marker birth file**", of the form **<run_ID>_birth.cdf<outtim>**, which contains the deposition data (position, energy, pitch...) for each marker deposited in the NUBEAM timeslice corresponding to <OUTTIM>.


Heavy use is made of the TRANSP plasma state file, which is produced if the TRANSP namelist setting "FI_OUTTIM" is set. 
 The file can then be found in the tar.gz files of the form <run_ID>_FI_TAR.GZ<FI_outtim>, for example: 134020D17_FI_TAR.GZ1
 Inside this tar.gz file is another tar file called 'tmp_tar.tar', inside which is the desired plasma state file '<run_ID>_ps_ts1_state.cdf'.
 
An example file path is: 
**D17/134020D17_FI_TAR.GZ1/tmp_tar/134020D17_ps_ts1_state.cdf**


Plenty of information about TRANSP CDF variables and namelist settings can be found at: https://users.euro-fusion.org/expert/transp/Webpages/files/TF.PLN (variables), https://users.euro-fusion.org/expert/transp/Webpages/index.html (general info),  https://w3.pppl.gov/~pshare/help/transp.htm (namelist settings & other info), https://w3.pppl.gov/~pshare/help/nubeam.htm (namelist settings & other info).


The **prepare_input.py** script calls the scripts in _inputhelpers_ to retrieve the Bfield, Efield, plasma, neutrals, wall, and marker data from the TRANSP files (2) & (3) listed above, and then outputs the ASCOT5 input file in the directory **runs**.

If the setting **check_input** is set to **True**, then the scripts in _inputchecks_ are called, which plot the profiles from A) the ASCOT5 file, B) the plasma state file, and C) the full CDF file. It is also checked numerically that marker data in the ASCOT5 file is identical to that in the marker birth file.

There is a small amount (maybe quantify this) of variation between the plasma state file and the full CDF profiles in the Efield case, due to the small differences between the 'PLFLX' and 'psipol' values in the full CDF and plasma state files respectively.

