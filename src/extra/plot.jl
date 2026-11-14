# Compatibility entry point for older notebooks. Graphical rendering is now
# shared with the optional PyPlot extension instead of a second implementation.
using BlueTangle
using PyPlot

plot_measurement(args...; kwargs...) = plotq(args...; backend=:pyplot, kwargs...)
