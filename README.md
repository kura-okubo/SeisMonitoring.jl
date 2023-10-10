# SeisMonitoring.jl

A Julia package of the tools for ambient noise seismology.

Documentation is available from the badge below:

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kura-okubo.github.io/SeisMonitoring.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kura-okubo.github.io/SeisMonitoring.jl/stable)

Other badges:

[![Documentation](https://github.com/kura-okubo/SeisMonitoring.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/kura-okubo/SeisMonitoring.jl/actions/workflows/documentation.yml)
[![Run tests](https://github.com/kura-okubo/SeisMonitoring.jl/actions/workflows/test.yml/badge.svg)](https://github.com/kura-okubo/SeisMonitoring.jl/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/kura-okubo/SeisMonitoring.jl/graph/badge.svg?token=iNq1WJH5bK)](https://codecov.io/gh/kura-okubo/SeisMonitoring.jl)
[![DOI](https://zenodo.org/badge/259752194.svg)](https://zenodo.org/badge/latestdoi/259752194)
[![Github All Releases](https://img.shields.io/github/downloads/kura-okubo/SeisMonitoring.jl/total.svg)]()

## Install SeisIO to Mac M1
Currently we have a problem when adding the `SeisIO.jl` using Mac M1 chip.
To avoid the error,
1. Use the Docker container following [**SeisMonitoring_Example**](https://github.com/kura-okubo/SeisMonitoring_Example).
2. Follow (https://github.com/jpjones76/SeisIO.jl/pull/94/commits/9a8b4510636e442e89f3d0a76f63abc56f1ab054).
We need to install the `SeisIO.jl` in develop directory:
```
]dev SeisIO 
```
Then, edit the `~/.julia/dev/SeisIO/Project.toml` such that
```
HDF5 = "0.12.3, 0.13, 0.14.2"
LightXML = "0.8.1, 0.9"
```
You can avoid the precompile error for those packages.

## Installation

Type the commands below in the Julia REPL:

```julia
using Pkg; Pkg.update();
Pkg.add(PackageSpec(name="SeisIO", version="1.2.1")); # Skip if you already installed as above
Pkg.add(PackageSpec(name="SeisNoise", version="0.5.3"));
Pkg.develop(url="https://github.com/kura-okubo/SeisDvv.jl");
Pkg.develop(url="https://github.com/kura-okubo/SeisMonitoring.jl");
```

## Tutorial
We created the notebook of the tutorial in the different github repository, [**SeisMonitoring_Example**](https://github.com/kura-okubo/SeisMonitoring_Example). You can find how to download the data, remove the transient signals, compute cross-correlations, stack the correlation functions and measure the dv/v.


You can access to the from the badge below:

<a href="https://nbviewer.org/github/kura-okubo/SeisMonitoring_Example/blob/main/code/run_seismonitoring.ipynb" target="_blank">
   <img align="left"
      src="https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.png"
      width="109" height="20">
</a>
<br><br>


or read the QR code:

<img src="/docs/src/assets/QRcode_seismonitoring_example.png" alt="QR" width="150"/>

See also the repoitory of [**SeisMonitoring_Paper**](https://github.com/kura-okubo/SeisMonitoring_Paper) for the post-processing using the dv/v over 20 years at Parkfield.

## Reference

