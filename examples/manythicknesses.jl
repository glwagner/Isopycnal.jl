using PyPlot, PyCall, NullableArrays

@pyimport numpy.ma as ma
PyObject(a::NullableArray) = pycall(ma.array, Any, a.values, mask=a.isnull)

function span(X)
  X1, X2 = extrema(X)
  abs(X2-X1)
end

if !isdefined(:loadvars) || loadvars 
  include("draketools.jl")

  lat, lon, topo, gam, X, Y, Z = getvars()
  x, y = X[:, :, 1], Y[:, :, 1]
  z = Z[1, 1, :]

  loadvars = false
end

Lx = span(X) 
Ly = span(Y) 
aspect = 62Lx/111Ly

ks = [10, 50, 75]
medgam = squeeze(median(gam, [1, 2]), (1, 2))

#gams = [27.6, 27.8, 28.0]
#dels = [0.05, 0.05, 0.02]

gams = 26.5:0.02:28.4
dels = 0.0001*ones(gams)

sigma = []
zgamma = []
r = []
for (i, gami) in enumerate(gams)
  zi = calcisopycnalsurface(gami, gam, z)
  sigmai = calcthickness(gami, gam, z; delta=dels[i])

  ri = (mean(sigmai[.!isnan.(sigmai)].^2)
    * mean(sigmai[.!isnan.(sigmai)].^(-2) ))

  push!(sigma, sigmai)
  push!(zgamma, zi)
  push!(r, ri)

  println(gami, " ", mean(sigmai[.!isnan.(sigmai)]))
  println("ri = ", ri)
end

close("all")
fig, ax = subplots()
plot(gams, r, ".")


#=
close("all")
fig, axs = subplots(nrows=3, figsize=(6, 8))
cbs = []

for (i, ax) in enumerate(axs)
  sca(ax)

  pcolormesh(x, y, PyObject(NullableArray(sigma[i], isnan.(sigma[i]))), 
    vmin=limz[i][1], vmax=limz[i][2], cmap="YlGnBu_r")

  ax[:set_aspect](aspect, adjustable="box")
  ax[:set_xticks]([290, 295, 300, 305])
  lbl = [
    @sprintf("\$ \\gamma = %.2f \$", gams[i]),
    @sprintf("\$ \\Delta \\gamma = %.2f \$", dels[i]),
    @sprintf("\$ \\langle z_{\\gamma} \\rangle = %.1f \$ m", 
      mean(zgamma[i][.!isnan.(zgamma[i])]))
  ]
  text(290.2, -56.0, join(lbl, "\\\\[0.8ex] \n"), fontsize=8)

  cb = colorbar(shrink=0.8)
  #cb[:set_ticks](limz[i][1]:400:limz[i][2])
  push!(cbs, cb)


  #cb[:set_ticks](ticks[i])
end

axs[1][:tick_params](labelbottom=false, bottom=false, labeltop=true, top=true)
axs[2][:tick_params](labelbottom=false, bottom=false)

cbs[1][:ax][:xaxis][:set_label_position]("top") 
cbs[1][:ax][:set_xlabel](L"\sigma \, \mathrm{(m)}", fontsize=12, labelpad=8.0)

tight_layout()

savefig("thicknessmap.png", dpi=240)
=#
