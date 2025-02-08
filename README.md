# Blockwise Direct Search (BDS)

BDS is a package for solving nonlinear optimization problems without using derivatives. The current version can handle unconstrained problems. 

## What is BDS?

BDS is a derivative-free package using blockwise direct-search methods. The current version is implemented in MATLAB, and it is being implemented in other programming languages.

See [Haitian LI's presentation](https://lht97.github.io/documents/DFOS2024.pdf) on BDS for more information.

## How to install BDS?

1. Clone this repository. You should then get a folder named `bds` containing this README file and the
[`setup.m`](https://github.com/blockwise-direct-search/bds/blob/main/setup.m) file.

2. In the command window of MATLAB, change your directory to the above-mentioned folder, and execute

```matlab
setup
```

If the above succeeds, then the package `bds` is installed and ready to use. Try `help bds` for more information.

We do not support MATLAB R2017a or earlier. If there exists any problems, please open an issue by
https://github.com/blockwise-direct-search/bds/issues.

## The coverage of unit test (offered by [Codecov](https://about.codecov.io/))

[![Codecov](https://img.shields.io/codecov/c/github/blockwise-direct-search/bds?style=for-the-badge&logo=codecov)](https://app.codecov.io/github/blockwise-direct-search/bds)

## Test of BDS.
The tests are **automated** by [GitHub Actions](https://docs.github.com/en/actions).
- [![Check Spelling](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml)
- [![Unit test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml)
- [![Coverage test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test_coverage.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test_coverage.yml)
- [![Stress test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml)
- [![Parallel test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/parallel_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/parallel_test.yml)
- [![Recursive test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/recursive_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/recursive_test.yml)
- [![Verify norma](https://github.com/zeroth-order-optimization/bds/actions/workflows/verify_norma.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/verify_norma.yml)

The following tests are implemented by [Optiprofiler](https://github.com/optiprofiler/optiprofiler).
- [Tests](https://github.com/zeroth-order-optimization/bds/actions) at [zeroth-order-optimization/bds](https://github.com/zeroth-order-optimization/bds)

    - [![Profile cbds and original cbds using optiprofiler, big](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_orig_cbds_big.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_orig_cbds_big.yml)
    - [![Profile cbds and original cbds using optiprofiler, small](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_orig_cbds_small.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_orig_cbds_small.yml)
    - [![Profile cbds and ds using optiprofiler, big](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_ds_big.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_ds_big.yml)
    - [![Profile cbds and ds using optiprofiler, small](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_ds_small.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_ds_small.yml)
    - [![Profile original cbds and ds using optiprofiler, big](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_orig_cbds_ds_big.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_orig_cbds_ds_big.yml)
    - [![Profile original cbds and ds using optiprofiler, small](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_ds_small.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_ds_small.yml)
    - [![Profile cbds using different stopping criteria, big](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_stopping_criteria_big.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_stopping_criteria_big.yml)
    - [![Profile cbds using different stopping criteria, small](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_stopping_criteria_small.yml/badge.svg)](https://github.com/zeroth-order-optimization/bds/actions/workflows/profile_cbds_stopping_criteria_small.yml)

- [Tests](https://github.com/0thopt/bds/actions) at [0thopt/bds](https://github.com/0thopt/bds)

    - [![Profile cbds and bfo using optiprofiler, big](https://github.com/0thopt/bds/actions/workflows/profile_cbds_bfo_big.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_cbds_bfo_big.yml)
    - [![Profile cbds and bfo using optiprofiler, small](https://github.com/0thopt/bds/actions/workflows/profile_cbds_bfo_small.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_cbds_bfo_small.yml)
    - [![Profile original cbds and bfo using optiprofiler, big](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_bfo_big.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_bfo_big.yml)   
    - [![Profile original cbds and bfo using optiprofiler, small](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_bfo_small.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_bfo_small.yml) 
    - [![Profile cbds and lam using optiprofiler, big](https://github.com/0thopt/bds/actions/workflows/profile_cbds_lam_big.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_cbds_lam_big.yml)
    - [![Profile cbds and lam using optiprofiler, small](https://github.com/0thopt/bds/actions/workflows/profile_cbds_lam_small.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_cbds_lam_small.yml)
    - [![Profile original cbds and lam using optiprofiler, big](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_lam_big.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_lam_big.yml)
    - [![Profile original cbds and lam using optiprofiler, small](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_lam_small.yml/badge.svg)](https://github.com/0thopt/bds/actions/workflows/profile_orig_cbds_lam_small.yml)

- [Tests](https://github.com/bladesopt/bds/actions) at [bladesopt/bds](https://github.com/bladesopt/bds)

    - [![Profile cbds and bfgs using optiprofiler, big](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_bfgs_big.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_bfgs_big.yml)
    - [![Profile cbds and bfgs using optiprofiler, small](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_bfgs_small.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_bfgs_small.yml)
    - [![Profile original cbds and bfgs using optiprofiler, big](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_bfgs_big.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_bfgs_big.yml)
    - [![Profile original cbds and bfgs using optiprofiler, small](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_bfgs_small.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_bfgs_small.yml)
    - [![Profile cbds and newuoa using optiprofiler, big](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_newuoa_big.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_newuoa_big.yml)
    - [![Profile cbds and newuoa using optiprofiler, small](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_newuoa_small.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_cbds_newuoa_small.yml)
    - [![Profile original cbds and newuoa using optiprofiler, big](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_newuoa_big.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_newuoa_big.yml)
    - [![Profile original cbds and newuoa using optiprofiler, small](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_newuoa_small.yml/badge.svg)](https://github.com/bladesopt/bds/actions/workflows/profile_orig_cbds_newuoa_small.yml)

- [Tests](https://github.com/derivative-free-optimization/bds/actions) at [derivative-free-optimization/bds](https://github.com/derivative-free-optimization/bds)
    - [![Profile cbds with rotated initial point, big](https://github.com/derivative-free-optimization/bds/actions/workflows/profile_cbds_rotated_initial_point_big.yml/badge.svg)](https://github.com/derivative-free-optimization/bds/actions/workflows/profile_cbds_rotated_initial_point_big.yml)
    - [![Profile cbds with rotated initial point, small](https://github.com/derivative-free-optimization/bds/actions/workflows/profile_cbds_rotated_initial_point_small.yml/badge.svg)](https://github.com/derivative-free-optimization/bds/actions/workflows/profile_cbds_rotated_initial_point_small.yml)

  
- [Tests](https://github.com/dfopt/bds/actions) at [dfopt/bds](https://github.com/dfopt/bds)

    - [![Profile cbds with cycling all blocks, big](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_all_big.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_all_big.yml)
    - [![Profile cbds with cycling all blocks, small](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_all_small.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_all_small.yml) 
    - [![Profile cbds with cycling all blocks using strategy 1, big](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_1_big.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_1_big.yml)
    - [![Profile cbds with cycling all blocks using strategy 1, small](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_1_small.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_1_small.yml)  
    - [![Profile cbds with cycling all blocks using strategy 2, big](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_2_big.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_2_big.yml)
    - [![Profile cbds with cycling all blocks using strategy 2, small](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_2_small.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_2_small.yml)   
    - [![Profile cbds with cycling all blocks using strategy 3, big](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_3_big.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_3_big.yml)
    - [![Profile cbds with cycling all blocks using strategy 3, small](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_3_small.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_3_small.yml)    
    - [![Profile cbds with cycling all blocks using strategy 4, big](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_4_big.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_4_big.yml)
    - [![Profile cbds with cycling all blocks using strategy 4, small](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_4_small.yml/badge.svg)](https://github.com/dfopt/bds/actions/workflows/profile_cbds_cycle_blocks_single_4_small.yml)   

- [Tests](https://github.com/gradient-free-opt/bds/actions) at [gradient-free-opt/bds](https://github.com/gradient-free-opt/bds)

    - [![Profile cbds and simplex using optiprofiler, big](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_simplex_big.yml/badge.svg)](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_simplex_big.yml)
    - [![Profile cbds and simplex using optiprofiler, small](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_simplex_small.yml/badge.svg)](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_simplex_small.yml)
    - [![Profile original cbds and simplex using optiprofiler, big](https://github.com/gradient-free-opt/bds/actions/workflows/profile_orig_cbds_simplex_big.yml/badge.svg)](https://github.com/gradient-free-opt/bds/actions/workflows/profile_orig_cbds_simplex_big.yml)
    - [![Profile original cbds and simplex using optiprofiler, small](https://github.com/gradient-free-opt/bds/actions/workflows/profile_orig_cbds_simplex_small.yml/badge.svg)](https://github.com/gradient-free-opt/bds/actions/workflows/profile_orig_cbds_simplex_small.yml)
    - [![Profile cbds and pds using optiprofiler, big](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_pds_big.yml/badge.svg)](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_pds_big.yml)
    - [![Profile cbds and pds using optiprofiler, small](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_pds_small.yml/badge.svg)](https://github.com/gradient-free-opt/bds/actions/workflows/profile_cbds_pds_small.yml)

- [Tests](https://github.com/libblades/bds/actions) at [libblades/bds](https://github.com/libblades/bds)

    - [![Profile cbds and nomad, small](https://github.com/libblades/bds/actions/workflows/profile_cbds_nomad_small.yml/badge.svg)](https://github.com/libblades/bds/actions/workflows/profile_cbds_nomad_small.yml)
    - [![Profile original cbds and nomad, small](https://github.com/libblades/bds/actions/workflows/profile_orig_cbds_nomad_small.yml/badge.svg)](https://github.com/libblades/bds/actions/workflows/profile_orig_cbds_nomad_small.yml)
    - [![Profile cbds and nomad, big](https://github.com/libblades/bds/actions/workflows/profile_cbds_nomad_big.yml/badge.svg)](https://github.com/libblades/bds/actions/workflows/profile_cbds_nomad_big.yml)
    - [![Profile original cbds and nomad, big](https://github.com/libblades/bds/actions/workflows/profile_orig_cbds_nomad_big.yml/badge.svg)](https://github.com/libblades/bds/actions/workflows/profile_orig_cbds_nomad_big.yml)      

- [Tests](https://github.com/opt-lab/bds/actions) at [opt-lab/bds](https://github.com/opt-lab/bds)

    - [![Profile cbds with randomized orthogonal matrix input, big](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_orthogonal_big.yml/badge.svg)](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_orthogonal_big.yml)
    - [![Profile cbds with randomized orthogonal matrix input, small](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_orthogonal_small.yml/badge.svg)](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_orthogonal_small.yml)
    - [![Profile cbds with randomized gaussian matrix input, big](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_gaussian_big.yml/badge.svg)](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_gaussian_big.yml)
    - [![Profile cbds with randomized gaussian matrix input, small](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_gaussian_small.yml/badge.svg)](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_randomized_gaussian_small.yml)
    - [![Profile cbds with permuted matrix input, big](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_permuted_big.yml/badge.svg)](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_permuted_big.yml)
    - [![Profile cbds with permuted matrix input, small](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_permuted_small.yml/badge.svg)](https://github.com/opt-lab/bds/actions/workflows/profile_cbds_permuted_small.yml)


- [Tests](https://github.com/optimlib/bds/actions) at [optimlib/bds](https://github.com/optimlib/bds)
 
    - [![Profile cbds with different number of blocks, big](https://github.com/optimlib/bds/actions/workflows/profile_cbds_block_number_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_cbds_block_number_big.yml)
    - [![Profile cbds with different number of blocks, small](https://github.com/optimlib/bds/actions/workflows/profile_cbds_block_number_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_cbds_block_number_small.yml)
    - [![Profile rbds with different number of num_selected_blocks, big](https://github.com/optimlib/bds/actions/workflows/profile_rbds_num_selected_blocks_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_rbds_num_selected_blocks_big.yml)
    - [![Profile rbds with different number of num_selected_blocks, small](https://github.com/optimlib/bds/actions/workflows/profile_rbds_num_selected_blocks_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_rbds_num_selected_blocks_small.yml)
    - [![Profile rbds with different replacement_delay, big](https://github.com/optimlib/bds/actions/workflows/profile_rbds_replacement_delay_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_rbds_replacement_delay_big.yml)
    - [![Profile rbds with different replacement_delay, small](https://github.com/optimlib/bds/actions/workflows/profile_rbds_replacement_delay_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_rbds_replacement_delay_small.yml)
    - [![Profile rbds with nsb n and rbds with nsb 1 and delay n-1, big](https://github.com/optimlib/bds/actions/workflows/profile_rbds_nsb_n_rbds_nsb_1_delay_n-1_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_rbds_nsb_n_rbds_nsb_1_delay_n-1_big.yml) 
    - [![Profile rbds with nsb n and rbds with nsb 1 and delay n-1, small](https://github.com/optimlib/bds/actions/workflows/profile_rbds_nsb_n_rbds_nsb_1_delay_n-1_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_rbds_nsb_n_rbds_nsb_1_delay_n-1_small.yml) 
    - [![Profile pbds and rbds with nsb n, big](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_n_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_n_big.yml) 
    - [![Profile pbds and rbds with nsb n, small](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_n_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_n_small.yml) 
    - [![Profile pbds and rbds with nsb 1 and delay n-1, big](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_1_delay_n-1_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_1_delay_n-1_big.yml)
    - [![Profile pbds and rbds with nsb 1 and delay n-1, small](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_1_delay_n-1_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_pbds_rbds_nsb_1_delay_n-1_small.yml)
    - [![Profile cbds and rbds with nsb n, big](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_n_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_n_big.yml) 
    - [![Profile cbds and rbds with nsb n, small](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_n_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_n_small.yml) 
    - [![Profile cbds and rbds with nsb 1 and delay n-1, big](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_1_delay_n-1_big.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_1_delay_n-1_big.yml)
    - [![Profile cbds and rbds with nsb 1 and delay n-1, small](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_1_delay_n-1_small.yml/badge.svg)](https://github.com/optimlib/bds/actions/workflows/profile_cbds_rbds_nsb_1_delay_n-1_small.yml)  
   

- [Tests](https://github.com/gradient-free-optimization/bds/actions) at [gradient-free-optimization/bds](https://github.com/gradient-free-optimization/bds)

    - [![Profile cbds and pads using optiprofiler, big](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pads_big.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pads_big.yml)
    - [![Profile cbds and pads using optiprofiler, small](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pads_small.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pads_small.yml) 
    - [![Profile cbds and scbds using optiprofiler, big](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_scbds_big.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_scbds_big.yml)
    - [![Profile cbds and scbds using optiprofiler, small](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_scbds_small.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_scbds_small.yml) 
    - [![Profile cbds and pbds using optiprofiler, big](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pbds_big.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pbds_big.yml)
    - [![Profile cbds and pbds using optiprofiler, small](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pbds_small.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_pbds_small.yml)  
    - [![Profile cbds and rbds using optiprofiler, big](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_rbds_big.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_rbds_big.yml)
    - [![Profile cbds and rbds using optiprofiler, small](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_rbds_small.yml/badge.svg)](https://github.com/gradient-free-optimization/bds/actions/workflows/profile_cbds_rbds_small.yml)   