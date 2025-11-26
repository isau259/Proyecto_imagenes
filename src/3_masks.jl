module Masks

export variable_density_mask

using Random

function variable_density_mask(Nx::Int, Ny::Int; accel::Real=4, center_fraction::Real=0.1, rng=Random.default_rng())
    mask = zeros(Float64, Nx, Ny)
    center_lines = Int(round(center_fraction * Ny))
    start_c = Int(round((Ny - center_lines) / 2))
    mask[:, start_c:(start_c + center_lines - 1)] .= 1.0
    remaining = [j for j in 1:Ny if !(start_c <= j <= start_c + center_lines - 1)]
    n_samples = Int(round(length(remaining) / accel))
    if n_samples > 0
        idx = collect(remaining)
        Random.shuffle!(rng, idx)
        sampled = idx[1:n_samples]
        mask[:, sampled] .= 1.0
    end
    return mask
end

end 
