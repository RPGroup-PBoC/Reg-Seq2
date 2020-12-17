# Twist order

Here will be a summary of all elements in the twist order. The files have to be `.csv` files, and the sequence column has to be labeled as `seq`. Add another
column labeled `primer_added`, which is either `True` if the primer pair is already attached to the sequence, or `False`, if the primers still have to be added.

| File      | Description | Reference Notebook | Creator |
| --------- | ----------- | ------------------ | ------- |
| lacI_sequences.csv | lacUV5+O1 mutants | `code/experimental_design/twist_order/lacI_titration/generate_sequences.ipynb`| Tom |
| lacUV5_tetOx_single_double_mutants.csv | lacUV5+tetOx single and double mutants | `code/experimental_design/twist_order/tetR_regulation/generate_sequences.ipynb`| Tom |
| natural_tet_promoters.csv | P_tetR1, P_tetR2, P_tetA mutated at 0.1 rate | `code/experimental_design/twist_order/tetR_regulation/generate_sequences.ipynb`| Tom |
| lacUV5_mutants.csv | various mutants with up to six mutations of *lacUV5* | `code/experimental_design/twist_order/lacUV5_mutants/generate_sequences.ipynb`| Tom |
| twist_sys_scrambles_2_16.csv | Systematic scrambles of 29 promoters with window sizes from 2-16 bp. | `code/experimental_design/twist_order/systematic_scrambles/twist_systematic_scrambles.ipynb`| Scott |
| twist_sys_scrambles_10.csv | Systematic scrambles of 29 promoters with window size 10 bp and 1 bp overlaps. | `code/experimental_design/twist_order/systematic_scrambles/twist_systematic_scrambles.ipynb`| Scott |
| twist_site_scrambles.csv | Scrambles of biocyc defined TF binding sites in 29 promoters. | `code/experimental_design/twist_order/site_scrambles/twist_site_scrambles.ipynb`| Scott |
| twist_orbit_TF_del_first_last_short.csv | Orbit deletion oligos for TFs. Genes shorter than 575bp with first/last amino acid strategy. | `code/experimental_design/twist_order/orbit/twist_orbit_TF_deletion.ipynb`| Scott |
| twist_orbit_TF_del_first_last_long.csv | Orbit deletion oligos for TFs. Genes longer than 575bp with first/last amino acid strategy. | `code/experimental_design/twist_order/orbit/twist_orbit_TF_deletion.ipynb`| Scott |
| twist_orbit_TF_del_avd_ovlp_short.csv | Orbit deletion oligos for TFs. Genes shorter than 575bp with avoid overlap strategy. | `code/experimental_design/twist_order/orbit/twist_orbit_TF_deletion.ipynb`| Scott |
| twist_orbit_TF_del_avd_ovlp_long.csv | Orbit deletion oligos for TFs. Genes longer than 575bp with avoid overlap strategy. | `code/experimental_design/twist_order/orbit/twist_orbit_TF_deletion.ipynb`| Scott |
