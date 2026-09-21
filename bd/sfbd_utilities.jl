
"""
Discretises real-valued x into grid. Returns the index bin x is assigned to.
"""
function bin(x, grid)
  i = 0
  while true
    i += 1
    if x <= grid[i]+step(grid)/2 || i >= length(grid)
      break
    end
  end
  return i
end

function is_in(trajectory_bin, searcharound_index, searchfor_index, window_size=0)
  for k in -window_size:window_size
    if searcharound_index+k >= 1 && searcharound_index+k <= length(trajectory_bin)
      if trajectory_bin[searcharound_index+k] == searchfor_index
        return true
      end
    end
  end
  return false
end

function all_sums!(result, body)
    for (i,el) in enumerate(body)
        result[i] = el
    end
    for (i1,el1) in enumerate(body), (i2,el2) in enumerate(body)
      if el1+el2 <= body[end] #&& i2 >= i1
        result[i1*length(body) + i2] = el1+el2
      end
    end
    # filter!(e->e!=-1,result)
    # unique!(result)
    # sort!(result)
    return result
end
function all_sums(body)
    result = zeros(Int, length(body)^2+length(body)) .- 1
    all_sums!(result,body)
    return result
end

function correlation(connected_nn, g2p_nn, M, time_coarsening_indices)
  correlation, disconnected = zeros(Float64, size(connected_nn)), zeros(Float64, size(connected_nn))
  time_coarsening_indices_sums = all_sums(time_coarsening_indices)
  for (i_tau_log, i_tau) in enumerate(time_coarsening_indices), (i_t1_log, i_t1) in enumerate(time_coarsening_indices)
    if i_tau+i_t1 <= time_coarsening_indices[end]
      i_tau_p_t1 = i_tau_log*length(time_coarsening_indices) + i_t1_log
      disconnected[i_t1_log, i_tau_log] = g2p_nn[1,i_tau_p_t1] * g2p_nn[2,i_t1_log] / M^2
      correlation[i_t1_log, i_tau_log] = connected_nn[i_t1_log, i_tau_log] / M - disconnected[i_t1_log, i_tau_log]
    end
  end
  return correlation,connected_nn/M,disconnected
end


"""
Returns a time grid that is a subset ('coarsened') of time_grid containing
coarsened_length points and discarding the first offset grid points.
Usually coarsened_length << length(time_grid).
The returned coarsened grid is part logarithmic and part linear, set by fraction.
"""
function coarsen_time(time_grid, coarsened_length::Int, offset; fraction=0.25)
    len_lin = Int(round(coarsened_length*fraction))
    time_coarsening_indices  = Int.(round.(10 .^((log10(offset):log10(length(time_grid))/(coarsened_length-len_lin+1):log10(length(time_grid)))) ))
    time_coarsening_indices2 = Int.(round.(range(offset, length(time_grid), length=len_lin)))
    append!(time_coarsening_indices, time_coarsening_indices2)
    sort!(time_coarsening_indices)
    unique!(time_coarsening_indices)
    return time_coarsening_indices
end
