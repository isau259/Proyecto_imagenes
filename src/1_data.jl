module Data

using Images
using FileIO
using ImageTransformations: imresize 

export gaussian_phantom, disk_phantom, shepp_logan_phantom, load_real_image

function gaussian_phantom(N::Int)
    x = range(-1, 1, length=N)
    img = [exp(-((xi^2 + yj^2) / 0.2^2)) for xi in x, yj in x]
    img ./= maximum(img)
    return img
end

function disk_phantom(N::Int)
    x = range(-1, 1, length=N)
    img = zeros(Float64, N, N)
    r0 = 0.25
    for i in 1:N, j in 1:N
        if x[i]^2 + x[j]^2 <= r0^2
            img[i, j] = 1.0
        end
    end
    return img
end

struct Ellipse
    A::Float64
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    phi::Float64
end

function shepp_logan_phantom(N::Int)
    ellipses = [
        Ellipse( 1.0,   0.69,  0.92,  0.0,    0.0,    0.0),
        Ellipse(-0.8,   0.6624,0.874, 0.0,    0.0,    0.0),
        Ellipse(-0.2,   0.11,  0.31,  0.22,   0.0,   -18*pi/180),
        Ellipse(-0.2,   0.16,  0.41, -0.22,   0.0,    18*pi/180),
        Ellipse( 0.1,   0.21,  0.25,  0.0,    0.35,   0.0),
        Ellipse( 0.1,   0.046, 0.046, 0.0,    0.1,    0.0),
        Ellipse( 0.1,   0.046, 0.046, 0.0,   -0.1,    0.0),
        Ellipse( 0.1,   0.046, 0.023,-0.08,  -0.605,  0.0),
        Ellipse( 0.1,   0.023, 0.023, 0.0,   -0.605,  0.0),
        Ellipse( 0.1,   0.023, 0.046, 0.06,  -0.605,  0.0)
    ]

    x = range(-1, 1, length=N)
    y = range(-1, 1, length=N)
    img = zeros(Float64, N, N)

    for e in ellipses
        ce = cos(e.phi)
        se = sin(e.phi)
        for i in 1:N, j in 1:N
            xr =  (x[j] - e.x0)*ce + (y[i] - e.y0)*se
            yr = -(x[j] - e.x0)*se + (y[i] - e.y0)*ce
            if (xr/e.a)^2 + (yr/e.b)^2 <= 1.0
                img[i, j] += e.A
            end
        end
    end

    img .-= minimum(img)
    img ./= maximum(img)
    return img
end

function load_real_image(path::AbstractString; N::Int=256)
    img = load(path)
    img_gray = Gray.(img)
    x = Float64.(img_gray)
    if size(x) != (N, N)
        x = imresize(x, (N, N))
    end

    mn, mx = extrema(x)
    if mx > mn
        x = (x .- mn) ./ (mx - mn)
    else
        x .= 0.0
    end
    return x
end

end 