Gallery of use-cases
====================

Gallery of use-cases of epygram, to be generated for epygram sphinx documentation.

Usage
-----

* Get input data: `make get_inputs`

Then the notebooks can be opened from a jupyter browser window, or they can be ran in batch:

* for the epygram documentation: `make nb4doc`

* as an extendedset of _unit_ tests: `make nb_check`
  In this case, they can also be ran from the `tests` directory of EPyGrAM

* WARNING: before committing, clear outputs of the notebooks not to save it in Git:
  `make clear_output` or even `make clean` to also clear checkpoints (ignored in Git anyway)

