EPyGrAM installation
====================

# 0. EVENTUAL PRE-REQUISITES
-------------------------
0.1 if you don't have packages 'footprints' and 'bronx' already installed,
	install them beforehand (ideally in the site-packages)
0.2 it is advised to uninstall older version before installing, i.e.:
     - remove install directory
     - remove epygram-related variables in .bash_profile/.bashrc
     - restart session 
0.3 de-tar new version EPyGrAM-x.x.x.tar into the directory of your choice,
    hereafter named $EPYGRAM_INSTALL_DIR


# 1. LINK SOURCES IN SITE-PACKAGES
------------------------------- 
1.1 go to the site-packages directory,
    usually ~/.local/lib/pythonX.Y/site_packages/ (X.Y depending on your python version)
    but also accessible in python through site.getusersitepackages()
1.2 link: ln -s $EPYGRAM_INSTALL_DIR EPyGrAM
1.3 copy path file: cp $EPYGRAM_INSTALL_DIR/_install/EPyGrAM.pth .
1.4 add [absolute path to site-packages]/EPyGrAM/apptools to $PATH
    (e.g. export PATH=$PATH:~/.local/lib/pythonX.Y/site_packages/EPyGrAM/apptools)


# 2. USER CUSTOMIZATION
--------------------
2.1 create directory $HOME/.epygram if necessary
2.2 (optional) user customization: copy
  - userconfig_empty.py into $HOME/.epygram/userconfig.py
  - sfxflddesc_mod.F90 into $HOME/.epygram/sfxflddesc_mod.F90
  - gribapi.def.0/ into $HOME/.epygram/gribapi.def.0/
