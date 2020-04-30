************
Introduction
************

SeisMonitoring.jl is developed to perform:

* Download data from earthquake data center
* Format data for SeisMonitoring.jl
* Pre-process raw data such as

    - Basic filtering

        + band-pass
        + detrend
        + taper
        + downsampling
    - Remove transient signals (earthquake events, tremors, artificial noise)
    - Frequency decomposition using **wavelet transform filter**
* Compute cross-correlation

    - Fast computing for all station pairs and auto-correlations.
    - High order cross_correlation (C3).
    - Allowing for normalizing options:

        + cross-coherence
        + deconvolution

* Stack noise cross-correlation functions

    - Flexible moving-window stacking duration


To achieve high-performance seismic data processing, SeisMonitoring.jl uses `SeisIO.jl <https://seisio.readthedocs.io/en/latest/?badge=latest>`_ and
`SeisNoise.jl <https://tclements.github.io/SeisNoise.jl/latest/>`_ .


Installation
============
| Since it is under developing, we recommend install it with ``dev`` command.
| From your terminal, type ``]`` and::

    dev https://github.com/kura-okubo/SeisMonitoring.jl.git


Overview
========
Among many available packages for ambient noise seismology, our challenges of SeisMonitoring.jl are following:

* Simplify process flow for variety of application

    - Flexible input parameters over a wide range of scale in time and space.
    - Transparency of scripts and functions makes it easy to extend its application.
    - User-friendly `Gtk <https://juliagraphics.github.io/Gtk.jl/latest/>`_ GUI accelerates the configuration of projects

* Optimizing computation

    - Parallelize all process with high scalability

* Apply advanced processing techniques to improve monitoring quality

    - Removing transient signals
    - Wavelet transform filtering
    - Variety of stacking option

* Review and plot on-going status with Julia command prompt

    - Review the progress of each process
    - Plot waveform of raw data and cross-correlation function
    - Plot quality of velocity change measurement
