EPyGrAM installation:
--------------------

A. automatic installation:

execute "install_epygram.py -b -e"

B. manual installation:

1.1 create directory $HOME/.epygram if necessary
1.2 copy user_Field_Dict_FA.csv and userconfig.py into $HOME/.epygram/
    (optional, for customization possibilities) 
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
