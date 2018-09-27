MATLAB Compiler

1. Prerequisites for Deployment 

. Verify the MATLAB Compiler Runtime (MCR) is installed and ensure you    
  have installed version 7.15.   

. If the MCR is not installed, run MCRInstaller, located in:

  <matlabroot>*\toolbox\compiler\deploy\win64\MCRInstaller.exe

For more information on the MCR Installer, see the MATLAB Compiler 
documentation.    


NOTE: You will need administrator right to run MCRInstaller. 


2. Files to Deploy and Package

Files to package for Standalone 
================================
-THAMDOF_v3_4.exe
-MCRInstaller.exe 
   -include when building component by clicking "Add MCR" link 
    in deploytool
-This readme file 

3. Definitions

For a complete list of product terminology, go to 
http://www.mathworks.com/help and select MATLAB Compiler.



* NOTE: <matlabroot> is the directory where MATLAB is installed on the target machine.


