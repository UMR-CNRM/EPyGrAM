FAQ --- Frequently Asked Questions
==================================

.. highlight:: python

-----------------------------------------------------------

+ **I have weird errors (MPI?) when I try to read a FA file...**

    => the :func:`epygram.init_env` should be called at the beginning of your
    script(s), to be sure the :mod:`arpifs4py` library do not use MPI, or uses
    OpenMP with an initialized environment.

-----------------------------------------------------------

+ **I have weird ecCodes errors...**

    => check environment variables ``$[ECCODES|GRIB]_[SAMPLES|DEFINITION]_PATH``
    and make sure they are consistent with your ecCodes install or unset them

-----------------------------------------------------------

+ **There is a field in my Surfex FA file, recognized as H2DField by** ``epygram``
  **whereas it is a meta-data field; how can I manage to read it ?**

    => copy ``$EPYGRAM_INSTALL_DIR/src/epygram/data/sfxflddesc_mod.F90`` into ``$HOME/.epygram/sfxflddesc_mod.F90``,
    and add your field in there in the manner of another meta-data field of the same kind
    And warn the ``epygram`` team about it, so that it will enter next version
    default ``sfxflddesc_mod.F90``.

-----------------------------------------------------------

+ **I need to overwrite the faFieldName.def ecCodes definition file for FA fieldnames conversion to GRIB2 parameters:**

    => copy ``$EPYGRAM_INSTALL_DIR/src/epygram/data/gribapi.def.??/`` as ``$HOME/.epygram/gribapi.def.99/``,
    then edit ``$HOME/.epygram/gribapi.def.99/grib2/localConcepts/lfpw/faFieldName.def``
    and set any field therein.

-----------------------------------------------------------

+ **My PC does not have enough memory to deal with these global spectral fields...**
  
    >>> Legendre spectral transforms need XX.XX MB
    >>> memory, while only YY.YY MB is available:
    >>> SWAPPING prevented !
  
    => either:
  
    - run your script on supercomputers using the alias
      ``s1batch='sbatch -N 1 -p normal256 --mem 250000 -t 00:30:00'``
      e.g.: ``s1batch myscript.py options -of -my --script``
    - use the ``fa_sp2gp.py`` tool on supercomputers (using the above alias) to convert your
      file to all-gridpoint then work in gridpoint space

-----------------------------------------------------------

+ **I have a brutal crash (e.g. SegFault) while working with some large size fields...**
   
    => before running Python, you need to raise the stack size limit: ``ulimit -s unlimited``

-----------------------------------------------------------

+ **How to hide LFI/FA messages ?**

    => either:
  
    - ``export FA4PY_MUTE=1``
    - edit ``$HOME/.epygram/userconfig.py`` and set ``FA_mute_FA4py = True``

-----------------------------------------------------------

+ **How to modify values of parameters or options of ``epygram`` *configuration* ?**

    => values of :mod:`epygram.config` can be modified in
    ``$HOME/.epygram/userconfig.py``

-----------------------------------------------------------

+ **How to modify the validity of a FA file ?**

    => e.g. to set the cumulativeduration to 0:
    ``fa_resource.modify_validity(basis=fa_resource.validity.get(), cumulativeduration=datetime.timedelta(0))``

-----------------------------------------------------------

**(to be continued...)**


