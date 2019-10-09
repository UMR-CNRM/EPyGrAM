EPyGrAM installation:
--------------------

0.1 (it is advised to uninstall older version before installing, i.e.:
     - remove install directory
     - remove epygram-related variables in .bash_profile/.bashrc
     - restart session
     ) 
0.2 (de-tar new version EPyGrAM-x.x.x.tar into install directory)


1. add the content of file 'epygram_environment' to your $HOME/.bash_profile,
   meanwhile modifying the installation directory.


2.1 create directory $HOME/.epygram if necessary
2.2 (optional) user customization: copy
  - userconfig_empty.py into $HOME/.epygram/userconfig.py
  - sfxflddesc_mod.F90 into $HOME/.epygram/sfxflddesc_mod.F90
  - gribapi.def.0/ into $HOME/.epygram/gribapi.def.0/


3. source $HOME/.bash_profile or open a new terminal
