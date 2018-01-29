""" 
Load variables from the Drake Passage simulation.
"""
function getvars(;
  filedir  = joinpath("..", "data"),
  filename = "neutraldensity_drakepassage.mat")

  filepath = joinpath(filedir, filename)
  println("Loading \"", filepath, "\"...")

  vars = matread(filepath)

   lat = vars["lat"]
   lon = vars["lon"]
  topo = vars["topo"]
     γ = vars["neutral_density"]
     X = vars["X"]
     Y = vars["Y"]
     Z = vars["Z"]

     x, y, z = X[:, :, 1], Y[:, :, 1], Z[1, 1, :]
  nx, ny, nz = size(X)
      Lx, Ly = span(X), span(Y) 

  lat, lon, topo, γ, X, Y, Z, x, y, z, nx, ny, nz, Lx, Ly
end


