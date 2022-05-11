
.PHONY: nb4doc nb_check get_inputs clear_output clean_ipynb_checkpoints

nb4doc:
	bash ./run_notebooks.sh in_place

nb_check:
	bash ./run_notebooks.sh check_errors

get_inputs:
	bash ./get_inputs.sh

clear_output:
	bash ./run_notebooks.sh clear_output

clean_ipynb_checkpoints:
	rm -rf */.ipynb_checkpoints

clean: clear_output clean_ipynb_checkpoints
