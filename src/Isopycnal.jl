__precompile__()

module Isopycnal

using MAT, PyCall, NullableArrays
@pyimport numpy.ma as ma
export nanmean, interpcolumnvalue, interpsurfacevalue, calcisopycnalheight,
       span, calcthickness!, calcthickness, getbcbz, get3dkappa,
       isopycnalanalysis

PyObject(a::NullableArray) = pycall(ma.array, Any, a.values, mask=a.isnull)

include("plottools.jl")

function span(X)
  X1, X2 = extrema(X)
  abs(X2-X1)
end


nanmean(f::Function, x, nanidx::BitArray) = try 
  mean(f, x[.!nanidx])
catch
  NaN
end

nanmean(x, nanidx::BitArray) = nanmean(x->x, x, nanidx)
nanmean(f::Function, x) = nanmean(f, x, isnan.(x))
nanmean(x) = nanmean(x->x, x)


""" 
Returns a 3D diffusivity distribution that decays exponentially as 
a height above the bottom with bottom value kap0 and decay scale d.
The diffusivity is set to NaN for Z values below the bottom.
"""
function get3dkappa(Z, topo; kap0=1e-3, d=500.0)
  κ₃ = zeros(Z)
  @. κ₃ = kap0 * exp(-(Z+topo)/d)
  @. κ₃[Z<-topo] = NaN

  κ₃
end


"""
Analyze the isopycnals in the list γᵢ.

This function returns
  Zt      : Z-grid with regions below topography turned into NaNs.
  κ₃      : 3D distribution of diffusivity used in the analysis.
  σᵢ      : List of 2D arrays of thickness on each isopycnal γᵢ.
  zᵢ      : List of 2D arrays of height each isopycnal γᵢ.
  κᵢ      : List of 2D arrays of diffusivity on each isopycnal γᵢ.
  σbar    : List of average thickness on each isopycnal.
  zbar    : List of average height on each isopycnal.
  κbar    : List of average diffusivity on each isopycnal.
  r_un    : List of "uncorrelated enhancement" calculated on each isopycnal.
  r_co    : List of "correlated enhancement" calculated on each isopycnal.
"""
function isopycnalanalysis(γᵢ, γ, Z, topo; delta=0.001)

  # Preliminary calculations
  nx, ny, nz = size(γ)

  κ₃ = get3dkappa(Z, topo)
  Zt = deepcopy(Z)
  @. Zt[ Z < -topo ] = NaN

  # Allocate
  σᵢ = [ zeros(nx, ny) for i=1:length(γᵢ) ]   # thicknesses
  zᵢ = [ zeros(nx, ny) for i=1:length(γᵢ) ]   # heights
  κᵢ = [ zeros(nx, ny) for i=1:length(γᵢ) ]   # diffusivity
  κbar = zeros(γᵢ)  # Isopycnal-averaged κ
  zbar = zeros(γᵢ)  # Isopycnal-averaged height
  σbar = zeros(γᵢ)  # Isopycnal-averaged thickness  
  r_un = zeros(γᵢ)  # "Uncorrelated" estimated strain-induced enhancement
  r_co = zeros(γᵢ)  # "Correlated" estimated strain-indued enhancement

  for i = 1:length(γᵢ)
    interpsurfacevalue!(zᵢ[i], γᵢ[i], γ, Zt)
    interpsurfacevalue!(κᵢ[i], γᵢ[i], γ, κ₃)

    calcthickness!(σᵢ[i], γᵢ[i], γ, Zt; delta=delta)

    σbar[i] = nanmean(σᵢ[i])
    κbar[i] = nanmean(κᵢ[i], isnan.(σᵢ[i]))
    zbar[i] = nanmean(zᵢ[i], isnan.(σᵢ[i]))

    r_un[i] = nanmean(x->x^2, σᵢ[i]) * nanmean(x->x^(-2), σᵢ[i])
    r_co[i] = nanmean(x->x^2, σᵢ[i]) * nanmean(κᵢ[i] ./ σᵢ[i].^2) / κbar[i]

    msg  = @sprintf("γ: %.2f, mean depth: % 8.2f, mean(σ): % 5.1f m, ",
      γᵢ[i], nanmean(zᵢ[i], isnan.(σᵢ[i])), nanmean(σᵢ[i]))

    msg *= @sprintf("r_un: % 6.2f, r_co: %.2f", 
      r_un[i], r_co[i])

    println(msg)
  end

  Zt, κ₃, σᵢ, zᵢ, κᵢ, σbar, zbar, κbar, r_un, r_co  
end


"""
Get the value of the vector v at which at which ρ = ρ★. vals and ρ must
be one-dimensional arrays (aka, in z) with the same dimension.
"""
function interpcolumnvalue(ρ★, ρ, v) 

  # Assume that first element is shallowest, so that z[1] > z[2] > z[3] > ...
  # For stable stratification, this implies ρ[1] < ρ[2] < ρ[3] < ...

  ssf = searchsortedfirst(ρ, ρ★)

  if ssf == 1
    return NaN
  elseif ssf > length(ρ)
    return NaN
  else  
    # Linear interpolation
    i₁, i₂ = ssf-1, ssf
    return v[i₁] + (v[i₂] - v[i₁]) / (ρ[i₂] - ρ[i₁]) * (ρ★ - ρ[i₁]) 
  end
end


"""
    interpsurfacevalue(ρ★, ρ, v)

Returns a 2D array v★ that gives v at the depth where ρ[i, j] = ρ★ for each 
horizontal coordinate i, j
"""
function interpsurfacevalue(ρ★, ρ, v)
  v★ = zeros(nx, ny)
  interpsurfacevalue!(v★, ρ★, ρ, v)
  v★
end

function interpsurfacevalue!(v★, ρ★, ρ, v)
  nx, ny, nz = size(ρ)
  for j = 1:ny, i=1:nx
    @views v★[i, j] = interpcolumnvalue(ρ★, ρ[i, j, :], v[i, j, :])
  end
  nothing
end


"""
Get a 2D array of heights z★[i, j] at which ρ[i, j] = ρ★
"""
function calcisopycnalheight(ρ★, ρ::Array{Float64,3}, z::Array{Float64,1})
  z★ = zeros(nx, ny)
  calcisopycnalheight!(z★, ρ★, ρ, z)
  z★
end

function calcisopycnalheight!(z★, ρ★, ρ::Array{Float64,3}, z::Array{Float64,1})
  nx, ny, nz = size(ρ)
  for j = 1:ny, i=1:nx
    @views z★[i, j] = interpcolumnvalue(ρ★, ρ[i, j, :], z)
  end
  nothing
end

"""
Get a 1D array of heights z★ at which ρ = ρ★.
"""
function calcisopycnalheight(ρ★, ρ::Array{Float64,2}, z::Array{Float64,1})
  nh, nz = size(ρ)
  z★ = zeros(nh)
  for i=1:nh
    @views z★[i] = interpcolumnvalue(ρ★, ρ[i, :], z)
  end
  z★
end



"""
Get the thickness of the isopycnal layer between ρ₁ and ρ₂.
"""
function calcthickness!(σ, ρ₁, ρ₂, ρ, Z::Array{Float64,3}; 
  z₁=nothing, z₂=nothing)
  if z₁ == nothing; z₁ = interpsurfacevalue(ρ₁, ρ, Z)
  else;             interpsurfacevalue!(z₁, ρ₁, ρ, Z)
  end

  if z₂ == nothing; z₂ = interpsurfacevalue(ρ₂, ρ, Z)
  else;             interpsurfacevalue!(z₂, ρ₂, ρ, Z)
  end
  @. σ = abs(z₁-z₂)
  nothing
end


"""
Get the thickness of the isopycnal layer between ρ₁ and ρ₂.
"""
function calcthickness!(σ, ρ₁, ρ₂, ρ, z::Array{Float64,1}; 
  z₁=nothing, z₂=nothing)
  if z₁ == nothing; z₁ = calcisopycnalheight(ρ₁, ρ, z)
  else;             calcisopycnalheight!(z₁, ρ₁, ρ, z)
  end

  if z₂ == nothing; z₂ = calcisopycnalheight(ρ₂, ρ, z)
  else;             calcisopycnalheight!(z₂, ρ₂, ρ, z)
  end

  @. σ = abs(z₁-z₂)

  nothing
end

function calcthickness(ρ₀, ρ, z::Array{Float64,1}; delta=0.1)
  ρ₁, ρ₂ = ρ₀-0.5*delta, ρ₀+0.5*delta
  σ = zeros(size(ρ)[1:end-1])
  calcthickness!(σ, ρ₁, ρ₂, ρ, z)
  σ
end

function calcthickness!(σ, ρ₀, ρ, z; delta=0.1, z₁=nothing, z₂=nothing)
  ρ₁, ρ₂ = ρ₀-0.5*delta, ρ₀+0.5*delta
  calcthickness!(σ, ρ₁, ρ₂, ρ, z; z₁=nothing, z₂=nothing)
end

function getbcbz(b::Array{Float64,3}, z)
  nx, ny, nz = size(b)
  bc = zeros(nx, ny, nz-1)
  bz = zeros(nx, ny, nz-1)

  # Ordering is z[1] < z[2] < ...
  @. @views bz = (b[:, :, 1:end-1] - b[:, :, 2:end]) / (z[1:end-1] - z[2:end])
  @. @views bc = 0.5*(b[:, :, 1:end-1] + b[:, :, 2:end])
  zc = 0.5*(z[2:end]+z[1:end-1]) # z value at cell center

  zc, bc, bz
end

end # module
