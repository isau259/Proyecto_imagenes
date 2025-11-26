
module Operators

export F_op, F_adj, M_op, M_adj

using FFTW

function F_op(x::AbstractArray)
    return fftshift(fft(fftshift(x))) / sqrt(length(x))
end

function F_adj(k::AbstractArray)
    return fftshift(ifft(fftshift(k))) * sqrt(length(k))
end

M_op(k, mask) = mask .* k
M_adj(k, mask) = mask .* k

end 
