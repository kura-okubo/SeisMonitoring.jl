# __precompile__()
# module SeisXcorrelation
#
# using SeisIO, SeisNoise, Dates, Statistics, Geodesy, FFTW, DSP, JLD2, Sockets, Distributed
# using Geodesics
#
# # compute cross-correlation
# include("seisxcorr.jl")
# # compute reference cross-correlation function and stacking (linear stack, selective stack, ...)
# include("seisstack.jl")
# # compute convergence rate of stacked cc
# include("stackconvergence.jl")
#
# end # module
