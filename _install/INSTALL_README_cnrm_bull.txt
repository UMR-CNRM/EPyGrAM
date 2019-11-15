EPyGrAM installation
====================

A. automatic installation (recommended)
-------------------------

execute "install_epygram.py -b -e --link_eccodes"
(install_epygram.py -h to choose options)
You may need to execute with both python2 and python3 if needed.


B. manual installation
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
