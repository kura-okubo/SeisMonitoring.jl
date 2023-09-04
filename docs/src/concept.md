# Development concept of SeisMonitoring.jl

## Overview
We have developed `SeisMonitoring.jl` to perform the TB scale ambient seismic noise processings. The aims of this package are the following:

1. Mass downloading the continuous seismic data, removing transient signals, cross-correlations, stacking and measurement of dv/v.
2. Options of the processing steps such as the band-pass filtering using the wavelet transform, selective stacking and the dv/v measurement of stretching and MWCS.
2. Paralleling the processes using **multi-node** in cluster.
3. Optimizing the use of memory and the frequency of I/O.

The kernels of `SeisMonitoring.jl` to handle the seismic waveforms are `SeisIO.jl` and `SeisNoise.jl`.


| Name  | Functions | Reference |
| --- | --- | --- |
| [SeisIO.jl](https://seisio.readthedocs.io/en/latest/?badge=latest) | Handle the seismic waveforms | Jones, J. P., Okubo, K., Clements, T., and Denolle, M. A. Seisio: a fast, efficient geophysical data architecture for the julia language. Seismol. Res. Lett., 91(4):2368--2377, 2020, doi:[10.1785/0220190295](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/91/4/2368/583741/SeisIO-A-Fast-Efficient-Geophysical-Data?redirectedFrom=fulltext).
| [SeisNoise.jl](https://tclements.github.io/SeisNoise.jl/latest/) | Tools for ambient seismic noise processings | Clements, T. and Denolle, M. A. Seisnoise.jl: ambient seismic noise cross correlation on the cpu and gpu in julia. Seismol. Res. Lett., 92(1):517--527, 2020, doi:[10.1785/0220200192](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/92/1/517/591402/SeisNoise-jl-Ambient-Seismic-Noise-Cross?redirectedFrom=fulltext).
|

We handle the seismic waveforms using the data structures such as [`SeisIO.SeisData`](https://seisio.readthedocs.io/en/latest/src/Help/tutorial.html) and [`SeisNoise.CorrData`](https://tclements.github.io/SeisNoise.jl/latest/types/#CorrData-Objects-for-ambient-noise-cross-correlations), which contains the meta data including the history of applied processes on the data, such as the filtering and tapering, and the waveforms.

The functions of `SeisMonitoring.jl` have been developed in the separated packages such as [`SeisDownload.jl`](https://github.com/kura-okubo/SeisDownload.jl) for the ease of maintenance of the packages and dependencies, which is suitable for the Julia. However, we decided to merge them into a single package mainly to allow for using all the packages in one line, `using SeisMonitoring`.

The dependencies are managed by Julia system, and are listed in Project.toml. We minimized the dependencies, in particular not using the python-related modules and MPI.jl.

## Process parallelization
Most of the processes are parallelized using Julia-native `Distributed.pmap` function as the tasks associated with the ambient seismic noise processing can be process parallelized, i.e. asynchronized parallelization by station pairs, days, and frequency bands without the communication across the workers. The `pmap` allows for the multi-node parallelization. You can perform it by

```julia
using Distributed
NP = 480 # for 56 cores * 10 nodes. Some margins of cores to increase the RAM per core.
addprocs(NP) # Add the workers
print("After adding the procs: NP=$(nprocs())");
@everywhere using SeisMonitoring # You need to redefine the packages in all the processors.
```

### Parallelization in computing cross-correlation
The most computationally expensive process is **to compute the cross-correlations.** It is not only due to the number of tasks; **the duplication of the file I/O** needs to be optimized.

!!! note "Example of cross-correlation in parallel"
    Let us do the correlations across the stations `STA1`, `STA2`, `STA3` after applying the FFT on the waveforms. Then,
    the tasks are `STA1-STA2`, `STA1-STA3`, `STA2-STA3`.
    - If you read the seismic data or FFT files from the disk in each task, you need to read all the `STA1`, `STA2` and `STA3` twice from the disk. This is redundant, which could cause the damage in the disk.
    - The frequency of file I/O is very important to secure the scratch system and disks. See the documentation by TACC [https://docs.tacc.utexas.edu/tutorials/managingio/] (https://docs.tacc.utexas.edu/tutorials/managingio/) to optimize it.
    - So, you want to store all the FFTs in the memory and distribute them to the tasks. However, if you send the FFT from master core hosting the tasks to the workers in the different node, it takes time to transfer the data. You can run this metric, but it is inefficient comparing to conduct the tasks within a node.
    - Therefore, we first parallelized the time window e.g. every 2 years, and submitted the jobs separately, and each job executes the cross-correlations of all possible station pairs parallelized with the number of cores in the node. This metric is optimized in the file I/O such that the file of seismic data is accessed only once, and is distributed via the memory within the node.
    - It is **not** practical to store the data in a single large file even in the HDF5-based `JLD2` format. The I/O from a single large file is not recommended when you conduct the process parallelization. Therefore, we separated the seismic data into a daily length, and also output the correlation functions into the small pieces of '.jld2' files, which is gathered in the post-processing.

We monitored the usage of CPUs and memory using [TACC Remora ](https://docs.tacc.utexas.edu/software/remora/). You can find the case study of scaling of parallelization [in SeisMonitoring\_Paper/Appx/Scaling\_Frontera](https://github.com/kura-okubo/SeisMonitoring_Paper/tree/develop/Appx/Scaling_Frontera).


## Compatible to the other application
The default file format and processing parameters are tuned to process the data for the High Resolution Seismic Network (HRSN) (doi:[10.7932/HRSN)](https://ncedc.org/bp_doi_metadata.html); however, the work flow can be generalized for any kind of the data set. You can also reuse the internal functions to filter the data.
We left the comments and some deprecated functions as the reference for future development.
 Please customize the functions under the [MIT License](https://github.com/kura-okubo/SeisMonitoring.jl/blob/dev_parallel_modified/LICENSE).
