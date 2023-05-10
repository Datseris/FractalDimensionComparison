# Estimating the fractal dimension: a comparative review and open source implementations

This code base is accompanying our paper with title as above. It is available on arXiv: https://arxiv.org/abs/2109.05937

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project.
It is created by George Datseris.

To (locally) reproduce this project, do the following:

1. Download this code base as is and extract the zip somewhere.
2. Open a Julia console and do:
    ```
    julia> using Pkg
    julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
    julia> Pkg.activate("path/to/this/project/folder")
    julia> Pkg.instantiate()
    ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

The most important folder of this repository is `main`. It contains the scripts that generate the figures of the main sections of the paper. The code has been set up in such a way that when running the scripts the necessary data are produced and saved. The most likely scenario for a researcher that wants to use this code base is to simply alter the variables in one of the scripts of `main`, which will produce comparison plots for e.g., different datasets or different metaparameters.

In conjunction to `main` there is the `src` folder that contains the actual function definitions for computing the relevant quantities. Note that functions that estimate fractal dimensions are not defined in `src`. They are formally published in the [FractalDimensions.jl](https://juliadynamics.github.io/FractalDimensions.jl/stable/) submodule of DynamicalSystems.jl.

Some more explanations on the code:

- The scripts have `data` variables. In `src/data_generation` they are transformed into datasets. That's where you can define new datasets.
- Experimental data are uploaded to `data/experimental` and are specified by name. You can add more text files there and qualify them using name as in the example scripts.
- `src/make_C_H.jl` defines infrastructure that automatically computes and saves correlation sums utilizing DrWatson. Similarly for `src/make_EVT.jl` but now for the extreme value theory local dimensions.
- `src/fractal_dim_fit.jl` contains files for fitting slopes in linear regions.

Feel free to also look into other folders (`srcipts`) which contain non-polished code used during the preparation of the manuscript.
