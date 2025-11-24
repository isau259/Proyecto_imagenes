using Images
using FFTW

include("1_data.jl")
include("2_operators.jl")
include("3_masks.jl")
include("4_ista.jl")
include("5_metrics.jl")

using .Data: gaussian_phantom, disk_phantom, shepp_logan_phantom, load_real_image
using .Operators: F_op, F_adj, M_op, M_adj
using .Masks: variable_density_mask
using .ISTA: ista_l1
using .Metrics: psnr, ssim_index

N = 256
@assert iseven(N) "N debe ser par para la transformada Haar 2D"
accel = 4.0             
center_fraction = 0.10  
λ = 0.02
maxiter = 200
α = 0.8
println("Generando imagen de referencia...")
use_real_image = false   # true = usar imagen real, false = usar phantom
if use_real_image
    img_path = "mri_slice.png"
    x_true = load_real_image(img_path; N=N)
else
    # Puedes elegir el phantom que quieras
    # x_true = gaussian_phantom(N)
    x_true = disk_phantom(N)
    # x_true = shepp_logan_phantom(N)
end

println("Calculando k-space completo...")
k_full = F_op(x_true)

println("Creando máscara de submuestreo (acceleración = $accel)...")
mask = variable_density_mask(N, N; accel=accel, center_fraction=center_fraction)
total_points    = N * N
acquired_points = count(!iszero, mask)    
frac_acquired   = acquired_points / total_points
R_eff           = total_points / acquired_points

println("Puntos totales en k-space: $total_points")
println("Puntos adquiridos: $acquired_points")
println("Fracción adquirida: $(round(frac_acquired * 100, digits=2)) %")
println("Factor de aceleración efectivo R_eff ≈ $(round(R_eff, digits=2))")
println("Aplicando máscara a k-space...")
k_sub = M_op(k_full, mask)
println("Reconstrucción zero-filled (IFFT)...")
t_zf = @elapsed begin
    x_zf = real(F_adj(k_sub))
end
println("Reconstrucción CS con ISTA-L1 (wavelets Haar)...")
t_ista = @elapsed begin
    x_cs = ista_l1(k_sub, mask; λ=λ, maxiter=maxiter, α=α, verbose=true)
end

function norm01(x)
    mn, mx = extrema(x)
    mx == mn && return zeros(size(x))
    return (x .- mn) ./ (mx - mn)
end

x_true_n = norm01(x_true)
x_zf_n   = norm01(x_zf)
x_cs_n   = norm01(x_cs)

println("Calculando métricas...")
psnr_zf   = psnr(x_zf_n, x_true_n)
psnr_cs   = psnr(x_cs_n, x_true_n)
ssim_zf_v = ssim_index(x_zf_n, x_true_n)
ssim_cs_v = ssim_index(x_cs_n, x_true_n)

println()
println("=== Resultados ===")
println("Zero-filled:  PSNR = $(round(psnr_zf, digits=2)) dB,  SSIM = $(round(ssim_zf_v, digits=4))")
println("CS (ISTA-L1): PSNR = $(round(psnr_cs, digits=2)) dB,  SSIM = $(round(ssim_cs_v, digits=4))")

println()
println("=== Tiempos de cómputo ===")
println("Zero-filled (IFFT): $(round(t_zf, digits=4)) s")
println("CS (ISTA + wavelets): $(round(t_ista, digits=4)) s")

println("Guardando imágenes PNG...")
save("x_true.png", colorview(Gray, x_true_n))
save("x_zf.png",   colorview(Gray, x_zf_n))
save("x_cs.png",   colorview(Gray, x_cs_n))

println("Listo. Imágenes guardadas en el directorio actual.")
