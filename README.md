# S&P_Global_Intern

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> S&P_Global_Intern

It is authored by Haoyun Li.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "S&P_Global_Intern"
```
which auto-activate the project and enable local path handling from DrWatson.

There are several different showcase experiments in the scripts folder. The first one is the most simple one with several layers, the second one is with baffles, and the third one is with faults. Conducting simulations on different topography is for comparaing the impact of the difference in geological structure on the CO2 injection plume behavior. The fourth one is with larger range of depth, but with a sealed reservoir. The later ones are all with larger range and more reasonable box size. The fifth one is with several layers. The sixth one is with an horizontal injection well. The sixth is for with baffles while the seventh is with the fault.

Then, we also try with the 3D simulation model. The 3D model is based on unstructured mesh.