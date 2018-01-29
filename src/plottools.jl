using PyPlot

function removespines(spines...)
  ax = gca()
  for spine in spines
    ax[:spines][spine][:set_visible](false)  
  end
  nothing
end

function tickparams(; kwargs...)
  ax = gca()
  ax[:tick_params](; kwargs...)
  nothing
end

function axisleft()
  ax = gca()
  ax[:tick_params](right=false, labelright=false, left=true, labelleft=true)
  ax[:yaxis][:set_label_position]("left")
  nothing
end

function axisright()
  ax = gca()
  ax[:tick_params](left=false, labelleft=false, right=true, labelright=true)
  ax[:yaxis][:set_label_position]("right")
  nothing
end


