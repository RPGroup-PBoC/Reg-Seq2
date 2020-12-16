# Twist order

Here will be a summary of all elements in the twist order. The files have to be `.csv` files, and the sequence column has to be labeled as `seq`. Add another
column labeled `primer_added`, which is either `True` if the primer pair is already attached to the sequence, or `False`, if the primers still have to be added.

| File      | Description | Reference Notebook| Creator|
| ----------- | ----------- | ----------- | ----------- |
| lacI_sequences.csv      | lacUV5+O1 mutants       | `code/experimental_design/twist_order/lacI_titration/generate_sequences.ipynb`| Tom |
| lacUV5_tetOx_single_double_mutants.csv      | lacUV5+tetOx single and double mutants       | `code/experimental_design/twist_order/tetR_regulation/generate_sequences.ipynb`| Tom |
| natural_tet_promoters.csv      | P_tetR1, P_tetR2, P_tetA mutated at 0.1 rate       | `code/experimental_design/twist_order/tetR_regulation/generate_sequences.ipynb`| Tom |
| lacUV5_mutants.csv      | various mutants with up to six mutations of *lacUV5*       | `code/experimental_design/twist_order/lacUV5_mutants/generate_sequences.ipynb`| Tom |
