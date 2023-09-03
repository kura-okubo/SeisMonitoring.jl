# SeisMonitoring.jl

A Julia package of the tools for ambient noise seismology.

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kura-okubo.github.io/SeisMonitoring.jl/dev)

## Installation

Type the commands below in the Julia REPL:

```julia
using Pkg; Pkg.update();
Pkg.add(PackageSpec(name="SeisIO", version="1.2.1"));
Pkg.add(PackageSpec(name="SeisNoise", version="0.5.3"));
Pkg.add(url="https://github.com/kura-okubo/SeisDvv.jl");
Pkg.add(url="https://github.com/kura-okubo/SeisMonitoring.jl");
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

