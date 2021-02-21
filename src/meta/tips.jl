const tips_cache = Dict{Integer, AbstractInterpolation}()

function tips(iso_id::Integer, temp::Number)
    if iso_id ∉ keys(tips_cache)        
        path = joinpath(module_path, "meta", "q", @sprintf("q%d.jld2", iso_id))
        if isfile(path) === false
            return isotopologue(iso_id)[1, :q_t0]
        end
        data = load(path, "data")          
        tips_cache[iso_id] = LinearInterpolation(data[:,1], data[:,2])
    end    
    tips_cache[iso_id](temp)
end