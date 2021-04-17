const tips_cache = Dict{Integer, Matrix{Float64}}()

@inline function tips(iso_id::Integer, temp::T)::Float64 where T <: AbstractFloat        
    data::Matrix{Float64} = get!(tips_cache, iso_id) do 
        path = joinpath(module_path, "meta", "q", @sprintf("q%d.jld2", iso_id))
        if isfile(path) === false
            d = [c_T_ref isotopologue(iso_id)[1, :q_t0]; (c_T_ref + 1.0) isotopologue(iso_id)[1, :q_t0]]            
        else        
            load(path, "data")
        end
    end
    ind_lo = max(1, min(size(data)[1] - 1, searchsortedfirst(@view(data[:,1]), temp) - 1))
    ind_hi = min(size(data)[1], ind_lo + 1)        

    # return linear interpolation
    return (data[ind_hi,2] - data[ind_lo,2]) / (data[ind_hi,1] - data[ind_lo,1]) * (temp - data[ind_lo,1]) + data[ind_lo,2]
end