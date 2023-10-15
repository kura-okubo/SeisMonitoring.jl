# SeisMonitoring.jl

This is the documentation of SeisMonitoring.jl

## Install SeisIO.jl to Mac M1
Currently we have a problem when adding the `SeisIO.jl` using Mac M1 chip.
To avoid the error, try one of the options below:

- Option 1. Use the forked version of SeisIO. Type the follwing command in the Julia REPL:
```
using Pkg; Pkg.add(url="https://github.com/kura-okubo/SeisIO.jl");
```

- Option 2. Use the Docker container following 
[**SeisMonitoring_Example**](https://github.com/kura-okubo/SeisMonitoring_Example). 

- Option 3.  Follow (https://github.com/jpjones76/SeisIO.jl/pull/94/commits/9a8b4510636e442e89f3d0a76f63abc56f1ab054).
We need to install the `SeisIO.jl` in develop directory:
```
]dev SeisIO 
```
Then, edit the `~/.julia/dev/SeisIO/Project.toml` such that
```
HDF5 = "0.12.3, 0.13, 0.14.2"
LightXML = "0.8.1, 0.9"
```
You can avoid the precompile error for those packages. We corrected these dependencies in the forked version used in the option 1 above.

## Installation

Type the commands below in the Julia REPL:

```julia
using Pkg; Pkg.update();
# Pkg.add(PackageSpec(name="SeisIO", version="1.2.1")); 
Pkg.add(url="https://github.com/kura-okubo/SeisIO.jl") # forked version used for the environment including Mac M1
Pkg.add(PackageSpec(name="SeisNoise", version="0.5.3"));
Pkg.develop(url="https://github.com/kura-okubo/SeisDvv.jl");
Pkg.develop(url="https://github.com/kura-okubo/SeisMonitoring.jl");
```

!!! note
    The metric above is a stable way to install the packages as the `SeisDvv.jl` and `SeisMonitoring.jl` are not registered to the [General registory](https://github.com/JuliaRegistries/General), and it may cause the issue when using `]add https://github.com/kura-okubo/SeisMonitoring.jl` in the Julia REPL.

!!! tip "Uninstall the package"
    Type `] rm SeisMonitoring` in the Julia REPL to remove the SeisMonitoring.

## Tutorial
We created the notebook of the tutorial in a different github repository, [**SeisMonitoring_Example**](https://github.com/kura-okubo/SeisMonitoring_Example). You can find how to download the data, remove the transient signals, compute cross-correlations, stack the correlation functions and measure the dv/v.


You can access to the from the badge below:

```@raw html
<a href="https://nbviewer.org/github/kura-okubo/SeisMonitoring_Example/blob/main/code/run_seismonitoring.ipynb" target="_blank">
   <img align="left"
      src="https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.png"
      width="109" height="20">
</a>
<br><br>
```

or read the QR code:

```@raw html
<img src="./assets/QRcode_seismonitoring_example.png" alt="QR" width="150"/>
```

See also the repository of [**SeisMonitoring_Paper**](https://github.com/kura-okubo/SeisMonitoring_Paper) for the post-processing using the dv/v over 20 years at Parkfield.

## Reference
