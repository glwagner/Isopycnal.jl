if !isdefined(:loadvars) || loadvars 
  include("draketools.jl")
  lat, lon, topo, gam, X, Y, Z = getvars()
  x, y = X[:, :, 1], Y[:, :, 1]
  z = Z[1, 1, :]
  Lx = span(X) 
  Ly = span(Y) 
  aspect = 62Lx/111Ly
  loadvars = false
end

#γs = [28.16, 28.18, 28.20]
γs = [27.93, 27.94, 27.95]
dels = 0.001*[1, 1, 1] #[0.01, 0.01, 0.01]
limz = [[30, 100] for i=1:3] #[ [0, 30], [10, 40], [10, 150] ]

σ, zγ, r = [], [], []
for (i, γi) in enumerate(γs)
  zγi = calcisopycnalheight(γi, gam, z)
  σi = calcthickness(γi, gam, z; delta=dels[i])
  ri = nanmean(x->x^2, σi)*nanmean(x->x^(-2), σi)

  push!(σ, σi)
  push!(zγ, zγi)
  push!(r, ri)

  @printf "γ: %.2f kg/m^3, mean(σ): % 6.1f m, r: %.2f\n" γi nanmean(σi) ri
end



close("all")
fig, axs = subplots(nrows=3, figsize=(6, 8))
cbs = []

for (i, ax) in enumerate(axs)
  sca(ax)

  pcolormesh(x, y, PyObject(NullableArray(σ[i], isnan.(σ[i]))), 
    vmin=limz[i][1], vmax=limz[i][2], cmap="YlGnBu_r")

  ax[:set_aspect](aspect, adjustable="box")
  ax[:set_xticks]([290, 295, 300, 305])
  lbl = [
    @sprintf("\$ \\gamma = %.2f \$", γs[i]),
    @sprintf("\$ \\Delta \\gamma = %.2f \$", dels[i]),
    @sprintf("\$ \\langle z_{\\gamma} \\rangle = %.1f \$ m", 
      mean(zγ[i][.!isnan.(zγ[i])]))
  ]
  text(290.2, -56.0, join(lbl, "\\\\[0.8ex] \n"), fontsize=8)

  cb = colorbar(shrink=0.8)
  push!(cbs, cb)
end

axs[1][:tick_params](labelbottom=false, bottom=false, labeltop=true, top=true)
axs[2][:tick_params](labelbottom=false, bottom=false)

cbs[1][:ax][:xaxis][:set_label_position]("top") 
cbs[1][:ax][:set_xlabel](L"\sigma \, \mathrm{(m)}", fontsize=12, labelpad=8.0)

tight_layout()

savefig("thicknessmap-2.png", dpi=480)
