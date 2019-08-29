EPyGrAM installation:
--------------------

A. automatic installation:

execute "install_epygram.py -b -e --link_eccodes"
(install_epygram.py -h to choose options)

B. manual installation:

1.1 create directory $HOME/.epygram if necessary
1.2 (optional) user customization: copy
  - userconfig_empty.py into $HOME/.epygram/userconfig.py
  - sfxflddesc_mod.F90 into $HOME/.epygram/sfxflddesc_mod.F90
  - gribapi.def.0/ into $HOME/.epygram/gribapi.def.0/
1.3 (symbolically) link /home/common/epygram/EPyGrAM with linkname "src"
    in directory $HOME/.epygram/

2.1 copy the content of file '[bullx|cnrm]_profile' to your environment,
    for instance in a $HOME/.epygram_profile, that you may load from your .bash_profile:
    if [ -f ~/.epygram_profile ]; then
      . ~/.epygram_profile
    fi
2.2 if you want to point to a more often updated epygram,
    change the 'src' link to EPyGrAM-next

3. source $HOME/.epygram_profile or open a new session, if necessary
