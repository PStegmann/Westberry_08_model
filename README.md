# Description
This Fortran code is an implementation of the depth-resolved phytoplankton primary production model by Westberry et al. 2008.

# Getting Started

## Requirements
This code has been tested with gfortran 8.2.0 and Python 2.7.
The `plot.py` script requires NumPy and Matplotlib.

## Compiling the Code
1. Set the environmental variable `FCC` to your Fortran compiler of choice.
2. Run `make` in the `./build` folder.

## Using the Code

### Preliminaries
- If you are using the code, please cite the corresponding paper by Stegmann et al., (2018) in any publication (journal paper, presentation, poster,...) where you are using the code's results.
- If you plan to develop or modify the code, please create a `feature` branch from `develop` or fork the repository.

### Running the Code
To run the program type `./wb08` in the `./build` folder. The program output will be dumped directly to the shell.
To plot the output of  `wb08`, pipe the data to `./build/output/output.txt` and use the provided `./scripts/plot.py`.

# LICENSE:
Please see the file `LICENSE.txt` for details.

# References
Westberry, T., M. J. Behrenfeld, D. A. Siegel, and E. Boss (2008): Carbon-based primary productivity modeling with vertically resolved photoacclimation. Global Biogeochm. Cycl. 22, GB2024, doi:10.1029/2007GB003078 .

Stegmann, P. G., B. Sun, J. Ding, P. Yang and X. Zhang (2008): Study of the Effects of Phytoplankton Morphology 1 and Vertical Profile on Lidar Attenuated Backscatter and Depolarization Ratio. J. Quant. Spec. & Rad. Trans. (Under Review).
