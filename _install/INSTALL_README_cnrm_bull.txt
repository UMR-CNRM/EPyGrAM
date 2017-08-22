EPyGrAM installation:
--------------------

1.1 copy the content of file 'epygram_environment_[bull|cnrm]' to your environment,
    for instance in a $HOME/.epygram_profile, that you may load from your .bash_profile:
    if [ -f ~/.epygram_profile ]; then
      . ~/.epygram_profile
    fi
1.2 if you want to point to a more often updated epygram,
    change therein the VERSION of epygram to '-next'

2.1 create directory $HOME/.epygram if necessary
2.2 copy user_Field_Dict_FA.csv and userconfig.py into $HOME/.epygram/
    (optional, for customization possibilities) 


3. source $HOME/.epygram_profile or open a new session, if necessary
