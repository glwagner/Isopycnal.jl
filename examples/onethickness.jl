using Isopycnal, PyPlot

if !isdefined(:lat)
  include("drakepassage.jl")
  lat, lon, topo, gam, X, Y, Z = getvars()
  x, y = X[:, :, 1], Y[:, :, 1]
  z = Z[1, 1, :]
end

gami = 27.0
zi = calcisopycnalsurface(gami, gam, z)

close("all")
fig, ax2 = subplots()
pcolormesh(x, y, zi)
pause(0.1)

