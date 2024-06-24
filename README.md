This repository contains code to reproduce all theoretical models for "Strong interactions and isospin symmetry breaking in a supermoiré lattice".

Installation Instructions:
1. Clone this repository
2. Install Julia 1.10 (e.g. using https://github.com/JuliaLang/juliaup.)
3. Install packages needed for this code:
- Call `julia --project=/path/to/the/folder/where/this/repository/was/cloned/`.
- Alternatively, navigate to the folder with the command line, then call `julia --project=.`
- Press "]" to enter the julia package manager.
- Type `instantiate` then press return to install required packages.
- Check installation by exiting the package manager with `delete`/`backspace`, then call `using HD0d`.
4. Install Jupyter Kernel following instructions at https://julialang.github.io/IJulia.jl/stable/manual/installation/
- A kernel with multiple threads is recommended for speed.
5. Run the notebook `Examples/Supermoire_model_Figures.ipynb` to generate all theoretical data.
- Install the following packages for the notebook to run: `HDF5, Interpolations, DataFrames, ColorSchemes`
- For plotting, install PyPlot following instructions at `https://github.com/JuliaPy/PyPlot.jl`.
- The notebook with generate several folders inside the `Examples` folder to store data.
- The code takes about 10 minutes to run on a laptop in "quick" mode, but 30-60 minutes to
produce the figures used in the paper at full resolution.

Much of this code was adapted from https://github.com/erezberg/Dirac_revivals_theory that
was used for the theory of https://www.nature.com/articles/s41586-020-2373-y
