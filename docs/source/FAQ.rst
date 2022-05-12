FAQ --- Frequently Asked Questions
==================================

.. highlight:: python

-----------------------------------------------------------

+ **ecCodes does not seem to be found by Python**

    => (re-)install eccodes support for Python : ``pip3 install --user eccodes``

    if that is not enough, reinstall epygram:

    - clean any reference to ``epygram`` in your environment
    - ``unset PYTHONPATH``
    - @ CNRM, reinstall epygram running ``$EPYGRAM_INSTALL_DIR/_install/install_epygram.py -e``

-----------------------------------------------------------

+ **I have rather esoterical errors when I try to read a FA file...**

    => the :func:`epygram.init_env` may be initialized at the beginning of your
    script(s), to be sure the :mod:`arpifs4py` library do not use MPI, or uses
    OpenMP with an initialized environment.

-----------------------------------------------------------

+ **There is a field in my Surfex FA file, recognized as MiscField by** ``epygram``
  **whereas it is H2D, or vice-versa; how can I manage to read it ?**

    => edit ``$HOME/.epygram/sfxflddesc_mod.F90``, and add your field in the
    manner of the examples given therein or in ``$EPYGRAM_INSTALL_DIR/epygram/data/sfxflddesc_mod.F90``.
    And warn the ``epygram`` team about it, so that it will enter next version
    default ``sfxflddesc_mod.F90``.

-----------------------------------------------------------

+ **I need to overwrite the faFieldName.def ecCodes definition file for FA fieldnames conversion to GRIB2 parameters:**

    => copy ``$HOME/.epygram/gribapi.def.0/`` as ``$HOME/.epygram/gribapi.def.99/``,
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
    - use the ``fa_sp2gp.py`` tool on Bull (using the above alias) to convert your
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

**(to be continued...)**


