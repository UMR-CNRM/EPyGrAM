
.PHONY: nb4doc nb_check get_inputs clear_output clear_ipynb_checkpoints clean

get_inputs:
	bash ./get_inputs.sh

nb4doc:
	bash ./run_notebooks.sh in_place

nb_check:
	bash ./run_notebooks.sh check_errors

clear_output:
	bash ./run_notebooks.sh clear_output

clear_ipynb_checkpoints:
	rm -rf */.ipynb_checkpoints

clean: clear_output clear_ipynb_checkpoints
