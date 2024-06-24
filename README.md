This repository contains code to reproduce all theoretical models for "Strong interactions and isospin symmetry breaking in a supermoiré lattice".

Installation Instructions:
1. Clone this repository
2. Install Julia 1.10 (e.g. using https://github.com/JuliaLang/juliaup.)
3. Install packages needed for this code:
- Call `julia --project=/path/to/the/folder/where/this/repository/was/cloned/`
- Alternatively, navigate to the folder with the command line, then call `julia --project=.`
- Press "]" to go to the package manager
- Type `instantiate"` then press return to install required packages
- Check installation by exiting the package manager with `delete`/`backspace`, then call `using HD0d`.
4. Install Jupyter Kernel following instructions at https://julialang.github.io/IJulia.jl/stable/manual/installation/
- A kernel with multiple threads is recommended for speed
5. Run the notebook `Examples/Supermoire_model_Figures.ipynb` to generate all theoretical data.
- Install the following packages for the notebook to run: `HDF5, Interpolations, DataFrames, ColorSchemes`
- For plotting, install PyPlot following https://github.com/JuliaPy/PyPlot.jl
- The notebook with generate several folders inside the `Examples` folder to store data
- It takes about 10-60 minutes to run on a laptop
