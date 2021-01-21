EPyGrAM installation
====================


A. Automatic installation (recommended)
---------------------------------------

* On supercomputers:
module load python

then

* (given an EPyGrAM repository $EPYGRAM_INSTALL_DIR), run:
python3 $EPYGRAM_INSTALL_DIR/_install/setup_epygram.py [-v][-b] [-h for help about arguments]
and cf. complementary information printed out.
You may need to execute with both python2 and python3 if needed.


B. Manual installation
----------------------

# 1. LINK SOURCES IN SITE-PACKAGES
------------------------------- 
1.1 go to the site-packages directory,
    usually ~/.local/lib/pythonX.Y/site_packages/ (X.Y depending on your python version)
    but also accessible in python through site.getusersitepackages()
1.2 link: ln -s /home/common/epygram/EPyGrAM EPyGrAM
1.3 copy path file: cp EPyGrAM/_install/EPyGrAM.pth .
1.4 add [absolute path to site-packages]/EPyGrAM/apptools to $PATH
    (e.g. export PATH=$PATH:~/.local/lib/pythonX.Y/site_packages/EPyGrAM/apptools)


# 2. USER CUSTOMIZATION
--------------------
2.1 create directory $HOME/.epygram if necessary
2.2 (optional) user customization:
	from the above site-packages EPyGrAM/_install/ directory, copy
  	- userconfig_empty.py into $HOME/.epygram/userconfig.py
  	- sfxflddesc_mod.F90 into $HOME/.epygram/sfxflddesc_mod.F90
  	- gribapi.def.0/ into $HOME/.epygram/gribapi.def.0/
