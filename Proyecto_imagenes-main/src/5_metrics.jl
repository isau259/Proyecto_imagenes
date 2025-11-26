module Metrics

export mse, psnr, ssim_index

using Statistics

function mse(x_rec, x_true)
    return mean((Float64.(x_rec) .- Float64.(x_true)).^2)
end

function psnr(x_rec, x_true; maxval::Float64=1.0)
    m = mse(x_rec, x_true)
    return 10 * log10(maxval^2 / m)
end


function ssim_index(x_rec, x_true; maxval::Float64=1.0)
    X = Float64.(x_rec)
    Y = Float64.(x_true)

    μx = mean(X)
    μy = mean(Y)

    σx2 = mean((X .- μx).^2)
    σy2 = mean((Y .- μy).^2)
    σxy = mean((X .- μx) .* (Y .- μy))

    C1 = (0.01 * maxval)^2
    C2 = (0.03 * maxval)^2

    num  = (2μx * μy + C1) * (2σxy + C2)
    den  = (μx^2 + μy^2 + C1) * (σx2 + σy2 + C2)

    return num / den
end

end
