module ISTA

export ista_l1

using LinearAlgebra
using ..Operators: F_op, F_adj, M_op, M_adj

soft_threshold(v, τ) = sign.(v) .* max.(abs.(v) .- τ, 0.0)

function haar1d_forward(x::AbstractVector)
    N = length(x)
    @assert iseven(N) "La longitud del vector debe ser par para Haar"
    hN = div(N, 2)
    a = similar(x, hN)
    d = similar(x, hN)
    invsqrt2 = 1 / sqrt(2)
    @inbounds for i in 1:hN
        x1 = x[2i - 1]
        x2 = x[2i]
        a[i] = (x1 + x2) * invsqrt2 
        d[i] = (x1 - x2) * invsqrt2 
    end
    return vcat(a, d)
end

function haar1d_inverse(c::AbstractVector)
    N = length(c)
    hN = div(N, 2)
    a = view(c, 1:hN)
    d = view(c, hN+1:N)
    x = similar(c)
    invsqrt2 = 1 / sqrt(2)
    @inbounds for i in 1:hN
        s = a[i]
        w = d[i]
        x[2i - 1] = (s + w) * invsqrt2
        x[2i]     = (s - w) * invsqrt2
    end
    return x
end

function haar2d_forward(img::AbstractMatrix)
    Nx, Ny = size(img)
    @assert iseven(Nx) && iseven(Ny) "Las dimensiones de la imagen deben ser pares para Haar 2D"
    tmp = similar(img)
    @inbounds for i in 1:Nx
        row = view(img, i, :)
        tmp[i, :] = haar1d_forward(row)
    end
    coeffs = similar(img)
    @inbounds for j in 1:Ny
        col = view(tmp, :, j)
        coeffs[:, j] = haar1d_forward(col)
    end
    return coeffs
end

function haar2d_inverse(coeffs::AbstractMatrix)
    Nx, Ny = size(coeffs)
    @assert iseven(Nx) && iseven(Ny) "Las dimensiones de los coeficientes deben ser pares para Haar 2D"
    tmp = similar(coeffs)
    @inbounds for j in 1:Ny
        col = view(coeffs, :, j)
        tmp[:, j] = haar1d_inverse(col)
    end
    img = similar(coeffs)
    @inbounds for i in 1:Nx
        row = view(tmp, i, :)
        img[i, :] = haar1d_inverse(row)
    end
    return img
end

W_op(x) = haar2d_forward(x)
W_adj(c) = haar2d_inverse(c)

function grad_f(x, mask, y)
    k = F_op(x)
    res = M_op(k, mask) .- y
    return real(F_adj(M_adj(res, mask)))
end

function ista_l1(y, mask; λ=0.01, maxiter::Int=100, α::Float64=0.8, x0=nothing, verbose::Bool=true)
    x = x0 === nothing ? real(F_adj(y)) : copy(x0)
    for it in 1:maxiter
        g = grad_f(x, mask, y)
        x_tmp = x .- α .* g
        c = W_op(x_tmp)
        c_thr = soft_threshold(c, α * λ)
        x_new = W_adj(c_thr)
        diff = norm(x_new - x) / max(norm(x), eps())
        x = x_new
        if verbose && (it % 10 == 0)
            @info "ISTA iter $it, diff = $diff"
        end
        if diff < 1e-4
            verbose && @info "Convergencia alcanzada en iteración $it"
            break
        end
    end
    return x
end

end 
