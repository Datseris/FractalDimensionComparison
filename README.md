# Estimating the fractal dimension: a comparative review and open source implementations

This code base is accompanying our paper with title as see above. It is available on arXiv:

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

The most important folder of this repository is `main`. It contains the scripts that generate the figures of the paper. The code has been set up in such a way that when running the scripts the necessary data are produced and saved.
The `src` folder contains the source code that heavily builds on the functions from DynamicalSystems.jl and DrWatson. The `data` variables inside the scripts can be changed accordingly to what is available in `src/data_generation.jl`.

Other folders contain scripts and code we used during the research phase of the project. They may not be as "clean" as the `main` folder.