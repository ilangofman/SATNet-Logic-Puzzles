# SATNet Model for Solving KenKen and Numbrix Puzzles

This repository contains a SATNet solver for KenKen and Numbrix puzzles. We tested the original SATNet model on the capabilities of more complicated puzzles and extended the model with ensemble learning through bagging and stacking.

## Usage

To run the models, please view the corresponding Jupyter notebook. 

For KenKen puzzles, the following configuration variables are important and vary between configurations:

- `GAME_SIZE`: the value of `n` in an `n x n` board.
- `max_cage_total`: the maximum total of a cage.
- `MAX_CELL_SIZE`: the maximum number of cells in a cage.
- `MASK_OUT_SINGLE`: boolean flag to mask out single cells from learning.

The following libraries are required: PyTorch, NumPy, os, and json.


## Acknowledgments

We would like to thank the original authors of the SATNet model for their contribution to the field of SAT solving with Neural Networks. The original work can be found:

Po-Wei Wang, Priya L. Donti, Bryan Wilder, and J. Zico Kolter. Satnet: Bridging deep learning and logical reasoning
using a differentiable satisfiability solver. In Proceedings of ICML’19, pages 6545–6554, 2019. URL https: [github.com/locuslab/SATNet.git.](https://github.com/locuslab/SATNet.git)
