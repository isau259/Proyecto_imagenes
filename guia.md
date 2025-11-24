# Módulo `1_data.jl` – Generación de imágenes de prueba para MRI

El módulo `Data` proporciona las **imágenes de referencia** (phantoms sintéticos o imágenes reales) sobre las que se simula la adquisición de MRI y se prueba la reconstrucción por *Compressed Sensing* (CS).

En términos del modelo matemático de CS, aquí se construye la imagen verdadera $(x_{\text{true}})$ que luego se usa para:

* generar el espacio (k) completo mediante la transformada de Fourier,
* submuestrear con una máscara,
* reconstruir con zero–filled o CS,
* y calcular métricas (PSNR, SSIM).

---

## 1. Encabezado del módulo

```julia
module Data

using Images
using FileIO
using ImageTransformations: imresize 

export gaussian_phantom, disk_phantom, shepp_logan_phantom, load_real_image
```

* Declara el módulo `Data`, que agrupa funciones relacionadas con **imágenes de prueba**.

* Importa paquetes:

  * `Images`: manejo de imágenes en Julia (tipos, escala de grises, etc.).
  * `FileIO`: permite usar `load(path)` para leer archivos de imagen (PNG, JPG, etc.).
  * `ImageTransformations: imresize`: función para redimensionar imágenes a un tamaño deseado.

* `export ...` indica qué funciones estarán visibles fuera del módulo:

  * `gaussian_phantom`
  * `disk_phantom`
  * `shepp_logan_phantom`
  * `load_real_image`


## 2. `gaussian_phantom(N)`: fantoma gaussiano suave

```julia
function gaussian_phantom(N::Int)
    x = range(-1, 1, length=N)
    img = [exp(-((xi^2 + yj^2) / 0.2^2)) for xi in x, yj in x]
    img ./= maximum(img)
    return img
end
```

* Define un grid 2D uniforme `x` en el intervalo $[-1,1]$ con `N` puntos.
* Construye una imagen `img` donde cada píxel vale:

  $\text{img}(x_i, y_j) = \exp\left(-\frac{x_i^2 + y_j^2}{0.2^2}\right)$

* Es decir, un **“bulto” gaussiano isotrópico** concentrado en el centro.
* Finalmente normaliza la imagen dividiendo por el valor máximo para que quede en el rango $[0,1]$.
* Es una forma simple de partir probando la implementación de Fourier y CS antes de usar phantoms más complejos.

---

## 3. `disk_phantom(N)`: disco con bordes definidos

```julia
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
```

* Crea una imagen de tamaño `N × N` llena de ceros (fondo negro).
* Define un radio `r0 = 0.25`.
* Recorre todos los píxeles y marca `1.0` (blanco) en aquellos que satisfacen:

  $x_i^2 + x_j^2 \le r_0^2$
* El resultado es un **disco blanco** centrado, con bordes bien definidos.

* Este fantoma permite visualizar bien:

  * artefactos de Gibbs,
  * cómo CS recupera contornos,
  * comparación entre zero–filled y CS en un caso sencillo.

---

## 4. `Ellipse` y `shepp_logan_phantom(N)`: fantoma anatómico

### Estructura `Ellipse`

```julia
struct Ellipse
    A::Float64
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    phi::Float64
end
```

Cada elipse está definida por:

* `A`: intensidad (puede ser positiva o negativa).
* `a`, `b`: semiejes mayor y menor.
* `x0`, `y0`: centro de la elipse en coordenadas normalizadas.
* `phi`: ángulo de rotación en radianes.

Esta estructura permite describir **componentes elípticas** dentro de la imagen, que luego se suman para formar el fantoma.

### Fantoma de Shepp–Logan

```julia
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
```

1. Define un conjunto de 10 elipses con parámetros clásicos del **phantom de Shepp–Logan** (versión modificada):

   * una elipse grande que representa la “cabeza”,
   * zonas oscuras internas,
   * estructuras que simulan “lóbulos”, “lentes”, etc.

2. Para cada elipse `e`:

   * Calcula coseno y seno del ángulo `phi`.
   * Recorre todos los píxeles `(i,j)`:

     * Traslada y rota el punto $(x_j, y_i)$ al sistema de la elipse:
       
       $\begin{aligned}
       xr &= (x_j - x_0)\cos\phi + (y_i - y_0)\sin\phi \
       yr &= -(x_j - x_0)\sin\phi + (y_i - y_0)\cos\phi
       \end{aligned}$
       
     * Comprueba si está dentro de la elipse:
       
       $\left(\frac{xr}{a}\right)^2 + \left(\frac{yr}{b}\right)^2 \le 1.$
       
     * Si está dentro, suma la intensidad `A` en ese píxel.

3. Después de sumar todas las elipses:

   * Resta el mínimo (`img .-= minimum(img)`) para que no haya valores negativos.
   * Divide por el máximo (`img ./= maximum(img)`) para normalizar a ([0,1]).

* El Shepp–Logan es un **phantom estándar en tomografía y MRI**.
* Representa una sección de cabeza con diferentes regiones de intensidad:

  * bordes suaves y duros,
  * zonas oscuras y claras,
  * estructuras internas complejas.
* Es ideal para:

  * probar reconstrucciones,
  * ver cómo CS maneja detalles anatómicos,
  * comparar métodos de forma controlada.

---
## 5. `load_real_image(path; N)`: imagen real normalizada

```julia
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
```

1. **Carga** una imagen desde disco:

   * `load(path)` usa `FileIO` y detecta el formato automáticamente.

2. **Convierte a escala de grises**:

   * `img_gray = Gray.(img)` transforma cualquier imagen RGB a un objeto `Gray`, es decir, un único canal de intensidad.

3. **Convierte a matriz `Float64`**:

   * `x = Float64.(img_gray)` obtiene una matriz 2D de intensidades.

4. **Redimensiona si hace falta**:

   * Si el tamaño original de la imagen no es `N × N`, se aplica `imresize(x, (N, N))` para que sea compatible con el resto del pipeline (FFT, máscara, etc.).

5. **Normaliza a ([0,1])**:

   * Calcula `mn, mx = extrema(x)`.
   * Si `mx > mn`, aplica:
     
     $x \leftarrow \frac{x - \text{mn}}{\text{mx} - \text{mn}}$
     
   * Si por alguna razón `mx == mn` (imagen constante), deja la imagen en ceros.

* Permite incorporar una **imagen real de MRI** al mismo flujo que los phantoms:

  * tamaño estándar (N×N),
  * rango de intensidades normalizado,
  * lista para pasar a Fourier y simular submuestreo.
* Así puedes demostrar que tu método de CS no solo funciona en phantoms ideales, sino también en datos que se parecen más a una imagen clínica.

---

# Módulo `2_operators.jl` – Modelo de adquisición en espacio *k*

El módulo `Operators` implementa los **operadores lineales básicos** del modelo de MRI que se usa en el proyecto:

$y = M F x$

donde:

* (x) es la imagen en el espacio de objeto (dominio espacial),
* (F) es la transformada de Fourier 2D,
* (M) es la máscara de submuestreo en el espacio (k),
* (y) son los datos adquiridos (muestras de (k)-space).

Este módulo define:

* cómo pasamos de imagen a espacio (k) y viceversa,
* cómo aplicamos la máscara que simula la adquisición acelerada.

Es la “parte física” del modelo en tu problema inverso de *Compressed Sensing*.

---

## 1. Encabezado del módulo

```julia
module Operators

export F_op, F_adj, M_op

using FFTW
```

* Declara el módulo `Operators`, que agrupa los operadores de Fourier y la máscara de muestreo.
* Importa `FFTW`, la biblioteca que implementa las FFT rápidas en Julia.
* Exporta tres funciones principales:

  * `F_op`: transformada de Fourier directa,
  * `F_adj`: transformada inversa (adjunta),
  * `M_op`: operador de máscara en espacio (k).

---

## 2. `F_op(x)`: transformada de Fourier 2D unitaria

```julia
function F_op(x::AbstractArray)
    return fftshift(fft(fftshift(x))) / sqrt(length(x))
end
```

* `fftshift(x)`: reordena la imagen para que el origen de frecuencias (0,0) quede en el centro del array.
* `fft(…)`: aplica la transformada de Fourier (en 1D o 2D, según la dimensión del array `x`).
* Segundo `fftshift`: vuelve a centrar el resultado, de modo que la componente DC quede en el centro del espacio (k).
* Divide por `sqrt(length(x))`:

  * `length(x)` es el número total de elementos (`N * N` para una imagen 2D),
  * este factor hace que la FFT sea prácticamente **unitaria**, es decir, preserva la energía (norma 2) entre imagen y espacio (k).

En términos matemáticos, esta función implementa el operador lineal:

$F: x \mapsto \hat{x}$

con la convención de centrado usada habitualmente en MRI (DC en el centro de la matriz de k-space).

* En MRI, los datos crudos se encuentran en el **espacio (k)**, no en el espacio de imagen.
* Para simular la adquisición, se necesita pasar de la imagen (x) al dominio de Fourier (F x).
* Esta función es el **modelo directo** (sin máscara) del escáner ideal.

---

## 3. `F_adj(k)`: transformada inversa (adjunta de F)

```julia
function F_adj(k::AbstractArray)
    return fftshift(ifft(fftshift(k))) * sqrt(length(k))
end
```

* Aplica `fftshift` al espacio (k) para centrar la componente de baja frecuencia.
* Ejecuta `ifft(…)`: transformada de Fourier inversa.
* Usa de nuevo `fftshift` para reordenar el resultado al sistema de coordenadas usual.
* Multiplica por `sqrt(length(k))` para que sea el adjunto de `F_op` bajo la norma estándar.

En términos de operadores, `F_adj` implementa (F^*), el **adjunto** de (F).
Como `F` se definió de forma unitaria (por el factor `1/√N`), `F_adj` coincide con la transformada de Fourier inversa apropiadamente escalada.

* Se usa para reconstruir imágenes desde espacio (k) (por ejemplo, en la reconstrucción zero–filled).
* En el algoritmo ISTA, el gradiente del término de fidelidad requiere calcular:

  $\nabla f(x) = F^* M^* (M F x - y)$

  donde `F_adj` aparece explícitamente.
* Es un bloque fundamental en cualquier método iterativo de reconstrucción: proyecta los residuos del espacio (k) de vuelta al espacio de imagen.

---

## 4. `M_op(k, mask) y M_adj(k, mask)`

```julia
M_op(k, mask) = mask .* k
M_adj(k, mask) = mask .* k
```
* `mask .* k` realiza un **producto elemento a elemento** entre:

  * `mask`: matriz binaria (0/1) que indica qué puntos del espacio (k) se adquieren,
  * `k`: matriz compleja con los valores de k-space.

El resultado:

* En las posiciones donde `mask == 1`, se mantiene el valor original de `k`.
* Donde `mask == 0`, el valor se anula (equivalente a “no medido” o “rellenado con cero”).

`M_adj` se define igual porque, al ser una máscara real y diagonal, su adjunto (M^*) es el mismo operador.

* Representa el proceso de **submuestreo en espacio (k)**:

  * En la práctica, un resonador acelerado por CS no adquiere todas las líneas de k-space.
  * `M_op` simula exactamente “qué líneas se adquirieron”.
* En el modelo ($y = M F x$):

  * `F_op(x)` produce el k-space completo.
  * `M_op(F_op(x), mask)` da las **muestras realmente disponibles** (`y`).
* En reconstrucción:

  * Zero–filled: haces `F_adj(M_op(k_full, mask))`.
  * CS: usas `M_op` y `M_adj` dentro del cálculo del gradiente para imponer fidelidad a los datos medidos.


---

# Módulo `3_masks.jl` – Máscaras de submuestreo en espacio *k*

El módulo `Masks` define cómo se **submuestrea el espacio (k)** en los experimentos. Esto es clave porque en *Compressed Sensing* (CS) la aceleración viene de **adquirir menos muestras**.

Este módulo implementa la máscara `variable_density_mask`: patrón pseudo–aleatorio con densidad mayor en el centro de (k), pensado para CS.

---

## 1. Encabezado del módulo

```julia
module Masks

export variable_density_mask

using Random
```

* Declara el módulo `Masks`.
* Importa `Random` para poder barajar índices de forma reproducible.
* Exporta `variable_density_mask` (la principal para CS).

---

## 2. `variable_density_mask(Nx, Ny; accel, center_fraction, rng)`

```julia
function variable_density_mask(
    Nx::Int, Ny::Int;
    accel::Real=4,
    center_fraction::Real=0.1,
    rng=Random.default_rng()
)
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
```

Genera una máscara de tamaño `Nx × Ny` (mismo tamaño que el k-space) que tiene:

* **Líneas centrales de $(k_y)$** totalmente muestreadas.
* Algunas columnas adicionales en la periferia elegidas **al azar**.

El resultado es una máscara con **mayor densidad en el centro** y muestreo más escaso y aleatorio en las altas frecuencias, lo que es típico en CS para MRI.

### 2.2. Detalle paso a paso

1. **Inicialización de la máscara**

   ```julia
   mask = zeros(Float64, Nx, Ny)
   ```

   * Crea una matriz `Nx × Ny` llena de ceros.
   * Un `0.0` significa “no se adquiere” esa muestra de k-space.
   * Luego se irán poniendo `1.0` en las posiciones adquiridas.

2. **Bloque central totalmente muestreado**

   ```julia
   center_lines = Int(round(center_fraction * Ny))
   start_c = Int(round((Ny - center_lines) / 2))
   mask[:, start_c:(start_c + center_lines - 1)] .= 1.0
   ```

   * `center_fraction` es la fracción de columnas (en dirección $(k_y)$) que se mantienen **todas**.

     * Ejemplo: `center_fraction = 0.1` → se muestrea el 10 % central de las líneas en $(k_y)$.
   * `center_lines` calcula cuántas columnas son.
   * `start_c` calcula dónde empieza ese bloque central.
   * `mask[:, start_c:…] .= 1.0` marca **todas las filas** de esas columnas como adquiridas.

   **Motivación física**:
   El centro de k-space corresponde a las **bajas frecuencias**, que contienen la mayoría de la energía y del contraste global de la imagen. En CS es muy común asegurar que esa región se mida completa para mantener robustez en la reconstrucción.

3. **Columnas restantes (periferia)**

   ```julia
   remaining = [j for j in 1:Ny if !(start_c <= j <= start_c + center_lines - 1)]
   ```

   * Construye una lista con los índices de columna que **no** pertenecen al bloque central.
   * Son las zonas de altas frecuencias donde el submuestreo puede ser más agresivo.

4. **Número de columnas a muestrear en la periferia**

   ```julia
   n_samples = Int(round(length(remaining) / accel))
   ```

   * `accel` es el parámetro de “aceleración nominal” (por ejemplo, 4).
   * Aquí se decide que, de todas las columnas restantes, solo se muestreará aproximadamente un `1/accel` de ellas.
   * Si `accel = 4` → se toman aproximadamente el 25 % de las columnas periféricas.

5. **Selección aleatoria de columnas periféricas**

   ```julia
   if n_samples > 0
       idx = collect(remaining)
       Random.shuffle!(rng, idx)
       sampled = idx[1:n_samples]
       mask[:, sampled] .= 1.0
   end
   ```

   * Convierte `remaining` a un vector indexable (`idx`).
   * Baraja los índices con `Random.shuffle!`, usando el generador `rng`.
   * Toma las primeras `n_samples` posiciones (`sampled`).
   * Marca esas columnas en la máscara como `1.0` (adquiridas).

   **Motivación en CS**:

   * La aleatoriedad introduce **incoherencia** en el patrón de muestreo, que es uno de los pilares de CS.
   * El hecho de que el centro esté densamente muestreado y la periferia más dispersa es un tipo de **muestreo de densidad variable**, muy común en MRI acelerada por CS.

6. **Retorno**

   ```julia
   return mask
   ```

   * Se devuelve una matriz `mask` con entradas 0/1.
   * En el pipeline del proyecto, esta máscara se usa con `M_op(k_full, mask)` para simular datos submuestreados.

---

# Módulo `4_ista.jl` – Reconstrucción CS con ISTA y wavelets Haar

El módulo `ISTA` implementa el **algoritmo de reconstrucción por *Compressed Sensing*** que usas en el proyecto.
En particular, resuelve (de forma aproximada) el problema:


$\min_x ; \frac{1}{2},| M F x - y |_2^2 ;+; \lambda , | W x |_1$

donde:

* (x): imagen a reconstruir,
* (F): transformada de Fourier (modelo de MRI),
* (M): máscara de submuestreo en espacio (k),
* (y): datos adquiridos (k-space submuestreado),
* (W): transformada wavelet Haar 2D,
* (\lambda): peso de regularización (controla sparsidad).

El método usado es **ISTA** (*Iterative Shrinkage-Thresholding Algorithm*), con regularización $(\ell_1)$ en la base wavelet Haar.

---

## 1. Encabezado del módulo

```julia
module ISTA

export ista_l1

using LinearAlgebra
using ..Operators: F_op, F_adj, M_op, M_adj
```

* Define el módulo `ISTA`.
* Exporta la función principal `ista_l1`, que será llamada desde `main_demo.jl`.
* Importa:

  * `LinearAlgebra` para normas, etc.
  * Desde `Operators`: `F_op`, `F_adj`, `M_op`, `M_adj`, es decir, el modelo de adquisición (M F x) y sus adjuntos.

Este módulo es el encargado de **resolver el problema inverso** dado el modelo directo ya definido en `Operators.jl`.

---

## 2. `soft_threshold(v, τ)`: operador proximal de la norma $(\ell_1)$

```julia
soft_threshold(v, τ) = sign.(v) .* max.(abs.(v) .- τ, 0.0)
```

### ¿Qué hace?

Implementa el operador de **umbralización suave** (“soft thresholding”):

$\text{soft}(v_i, \tau) =
\begin{cases}
v_i - \tau,\text{sign}(v_i), & |v_i| > \tau \
0, & |v_i| \le \tau
\end{cases}$

Se aplica elemento a elemento sobre el vector/matriz `v` con un umbral `τ`.

### ¿Por qué es importante?

* Es el **operador proximal** de la norma $(\ell_1)$.

* En ISTA, esta operación aparece en el paso:
  
  $x^{(k+1)} = \text{prox}_{\alpha\lambda |\cdot|_1}\left( x^{(k)} - \alpha \nabla f(x^{(k)}) \right)$

* En tu implementación, no se aplica directamente sobre la imagen, sino sobre los **coeficientes wavelet** (W x), lo que promueve sparsidad en esa base.

---

## 3. Wavelet Haar 1D: `haar1d_forward` y `haar1d_inverse`

### `haar1d_forward(x)`

```julia
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
```

* Toma un vector `x` de longitud par.
* Lo agrupa en pares `(x1, x2)`.
* Calcula:

  * `a[i] = (x1 + x2) / √2`  → coeficientes de **aproximación** (baja frecuencia).
  * `d[i] = (x1 - x2) / √2`  → coeficientes de **detalle** (alta frecuencia).
* Devuelve un vector concatenado `[a; d]`.

Es la transformada Haar 1D en una **escala** (sin multirresolución por niveles).

### `haar1d_inverse(c)`

```julia
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
```

* Separa el vector de coeficientes `c` en `a` (aproximación) y `d` (detalle).

* Reconstruye el vector original:

  $x_{2i-1} = \frac{s + w}{\sqrt{2}}, \quad
  x_{2i}   = \frac{s - w}{\sqrt{2}}.$

* Es la inversa exacta de `haar1d_forward`.

---

## 4. Wavelet Haar 2D: `haar2d_forward` y `haar2d_inverse`

### `haar2d_forward(img)`

```julia
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
```

* Aplica Haar 1D por **filas** y luego por **columnas**:

  1. Para cada fila `i`, reemplaza la fila por su transformada Haar 1D.
  2. Luego, para cada columna `j`, aplica Haar 1D sobre esa columna.
* El resultado `coeffs` es la transformada wavelet Haar 2D de `img`.

Se trata de una transformada separable:
Haar 2D = Haar 1D en x seguido de Haar 1D en y.

### `haar2d_inverse(coeffs)`

```julia
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
```

* Aplica primero la **inversa Haar 1D por columnas**, luego por **filas**.
* Reconstruye la imagen original a partir de los coeficientes wavelet.

---

## 5. Operadores wavelet: `W_op` y `W_adj`

```julia
W_op(x) = haar2d_forward(x)
W_adj(c) = haar2d_inverse(c)
```

* `W_op`: transformada Haar 2D → pasa de imagen en el dominio espacial a coeficientes wavelet.
* `W_adj`: transformada inversa Haar 2D → pasa de coeficientes wavelet a imagen.

En tu formulación de CS:

$| W x |_1$

es precisamente la norma $(\ell_1)$ de los coeficientes Haar de la imagen $(x)$, y `W_op` y `W_adj` son los operadores que permiten aplicar ISTA en ese espacio.

---

## 6. `grad_f(x, mask, y)`: gradiente del término de fidelidad

```julia
function grad_f(x, mask, y)
    k = F_op(x)
    res = M_op(k, mask) .- y
    return real(F_adj(M_adj(res, mask)))
end
```

### ¿Qué hace?

Calcula el gradiente de:

$f(x) = \frac{1}{2} | M F x - y |_2^2.$

Paso a paso:

1. `k = F_op(x)`

   * Calcula el k-space completo de la imagen actual.

2. `res = M_op(k, mask) .- y`

   * Aplica la máscara: `M_op(k, mask)` da el k-space submuestreado simulado.
   * Resta los datos medidos `y` → residual en espacio (k).

3. `M_adj(res, mask)`

   * Aplica el adjunto de la máscara (igual que la máscara en este caso).

4. `F_adj(...)`

   * Lleva el residual de vuelta al espacio de imagen con el adjunto de Fourier.

5. `real(...)`

   * Se queda con la parte real (las imágenes que usas son reales).

En términos matemáticos:

$\nabla f(x) = F^* M^* (M F x - y).$

---

## 7. `ista_l1(y, mask; ...)`: algoritmo ISTA con (\ell_1) en Haar

```julia
function ista_l1(
    y, mask;
    λ=0.01,
    maxiter::Int=100,
    α::Float64=0.8,
    x0=nothing,
    verbose::Bool=true
)
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
```

### 7.1. Parámetros

* `y`: datos medidos en espacio (k), ya submuestreados (`M_op(k_full, mask)`).
* `mask`: máscara de submuestreo (variable density).
* `λ`: parámetro de regularización (peso de la sparsidad).
* `maxiter`: número máximo de iteraciones.
* `α`: paso de gradiente (step size).
* `x0`: inicialización. Si es `nothing`, se usa la reconstrucción zero–filled (`F_adj(y)`).
* `verbose`: controla si se muestran mensajes de progreso.

### 7.2. Inicialización

```julia
x = x0 === nothing ? real(F_adj(y)) : copy(x0)
```

* Si no se pasa una imagen inicial, se usa la **IFFT del k-space submuestreado** → reconstrucción zero–filled como punto de partida.
* Si se da `x0`, se copia.

### 7.3. Iteración ISTA

En cada iteración:

1. **Cálculo del gradiente de fidelidad**

   ```julia
   g = grad_f(x, mask, y)
   ```

   Esto computa $(\nabla f(x))$ para el término de datos.

2. **Paso de gradiente**

   ```julia
   x_tmp = x .- α .* g
   ```

   Implementa:

   $x^{(k + 1/2)} = x^{(k)} - \alpha \nabla f(x^{(k)}).$

3. **Transformada wavelet + umbralización suave**

   ```julia
   c = W_op(x_tmp)
   c_thr = soft_threshold(c, α * λ)
   ```

   * `W_op(x_tmp)`: pasa a coeficientes wavelet.
   * `soft_threshold(c, α * λ)`: aplica la proximidad de la norma (\ell_1) en el espacio wavelet, con umbral `α * λ`.

   Matemáticamente, esto corresponde a:

   $c^{(k+1)} = \text{soft}(W x^{(k+1/2)}, \alpha \lambda).$

4. **Transformada inversa wavelet**

   ```julia
   x_new = W_adj(c_thr)
   ```

   * Reconstruye la imagen desde los coeficientes umbralizados.
   * Este es el nuevo `x` en el dominio espacial.

5. **Criterio de convergencia**

   ```julia
   diff = norm(x_new - x) / max(norm(x), eps())
   x = x_new

   if verbose && (it % 10 == 0)
       @info "ISTA iter $it, diff = $diff"
   end
   if diff < 1e-4
       verbose && @info "Convergencia alcanzada en iteración $it"
       break
   end
   ```

   * Calcula un cambio relativo entre iteraciones.
   * Si `diff < 1e-4`, se considera que el algoritmo **ha convergido** y se detiene.
   * Cada 10 iteraciones, si `verbose` es `true`, imprime el valor de `diff` para monitorear.

### 7.4. Resultado

```julia
return x
```

Devuelve la imagen reconstruida, que es la solución aproximada del problema:

$\min_x ; \tfrac12 |M F x - y|_2^2 + \lambda |W x|_1.$

---
Te armo el `metrics.md` igual que los otros, explicando qué hace y por qué.

---

# Módulo `Metrics.jl` – Métricas de calidad de imagen (MSE, PSNR, SSIM)

El módulo `Metrics` agrupa las funciones que permiten **cuantificar** qué tan buena es una reconstrucción de imagen en comparación con una imagen de referencia (x_{\text{true}}).

Las tres métricas implementadas son:

* `mse`: *Mean Squared Error* (error cuadrático medio),
* `psnr`: *Peak Signal-to-Noise Ratio*,
* `ssim_index`: *Structural Similarity Index Measure* (SSIM) global.

Estas métricas se usan en tu proyecto para comparar:

* reconstrucción zero–filled vs.
* reconstrucción por CS (ISTA + wavelets),

y así mostrar numéricamente que CS **mejora la calidad de imagen** usando el mismo conjunto reducido de datos de k-space.

---

## 1. Encabezado del módulo

```julia
module Metrics

export mse, psnr, ssim_index

using Statistics
```

* Declara el módulo `Metrics`.
* Exporta las funciones `mse`, `psnr` y `ssim_index`, que se usarán desde `main_demo.jl`.
* Importa `Statistics` para usar funciones como `mean`.

---

## 2. `mse(x_rec, x_true)`: error cuadrático medio

```julia
function mse(x_rec, x_true)
    return mean((Float64.(x_rec) .- Float64.(x_true)).^2)
end
```

### ¿Qué hace?

* Convierte ambas imágenes a `Float64` (por si vienen en otros tipos).
* Calcula la diferencia punto a punto: `x_rec - x_true`.
* Eleva al cuadrado cada diferencia.
* Toma el promedio de todos los píxeles (`mean`).

Formalmente:

$\text{MSE}(x_{\text{rec}}, x_{\text{true}}) =
\frac{1}{N} \sum_{i=1}^N \big(x_{\text{rec}, i} - x_{\text{true}, i}\big)^2,$

donde (N) es el número total de píxeles.

### ¿Por qué es útil?

* Es una medida básica de error global.
* Cuanto **más pequeño** es el MSE, mejor se aproxima la reconstrucción a la imagen verdadera.
* Es la métrica subyacente que se usa en la definición estándar de PSNR.

---

## 3. `psnr(x_rec, x_true; maxval=1.0)`: relación pico–ruido

```julia
function psnr(x_rec, x_true; maxval::Float64=1.0)
    m = mse(x_rec, x_true)
    return 10 * log10(maxval^2 / m)
end
```

### ¿Qué hace?

* Llama a `mse(x_rec, x_true)` para obtener el error cuadrático medio.
* Aplica la fórmula clásica de PSNR:

$\text{PSNR} = 10 \log_{10}\left( \frac{\text{maxval}^2}{\text{MSE}} \right),$

donde:

* `maxval` es el valor máximo posible en la imagen (por defecto 1.0, porque normalizas las imágenes a ([0,1])),
* `MSE` es el error cuadrático medio entre reconstrucción y referencia.

### Interpretación

* **Cuanto mayor es el PSNR**, mejor la reconstrucción (menos error).
* Se mide en decibelios (dB).
* Es muy utilizada en procesamiento de imágenes para comparar métodos de compresión, filtrado y reconstrucción.

En tu proyecto:

* Comparar PSNR de zero–filled vs. CS permite mostrar numéricamente cuánto mejora la reconstrucción con *Compressed Sensing* para una misma máscara de submuestreo.

---

## 4. `ssim_index(x_rec, x_true; maxval=1.0)`: índice de similitud estructural

```julia
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
```

### 4.1. ¿Qué es SSIM?

SSIM (*Structural Similarity Index Measure*) intenta medir la similitud entre dos imágenes desde el punto de vista de:

* **luminancia** (brillo medio),
* **contraste** (varianza),
* **estructura** (correlación).

Su valor:

* es típicamente $(\in [0, 1])$ para imágenes [0,1],
* **1** significa imágenes prácticamente idénticas,
* mientras más cercano a 1, mejor.

Aquí implementas una versión **global** (no por ventanas), es decir, usas medias y varianzas globales de toda la imagen.

### 4.2. Detalle paso a paso

1. **Conversión a `Float64`**

   ```julia
   X = Float64.(x_rec)
   Y = Float64.(x_true)
   ```

   Asegura que ambas imágenes estén en el mismo tipo numérico.

2. **Cálculo de medias**

   ```julia
   μx = mean(X)
   μy = mean(Y)
   ```

   * `μx`: luminancia media de la reconstrucción.
   * `μy`: luminancia media de la referencia.

3. **Varianzas y covarianza**

   ```julia
   σx2 = mean((X .- μx).^2)
   σy2 = mean((Y .- μy).^2)
   σxy = mean((X .- μx) .* (Y .- μy))
   ```

   * `σx2`: varianza de `X` (contraste de la reconstrucción).
   * `σy2`: varianza de `Y` (contraste de la referencia).
   * `σxy`: covarianza conjunta (estructura compartida entre ambas imágenes).

4. **Constantes de estabilización**

   ```julia
   C1 = (0.01 * maxval)^2
   C2 = (0.03 * maxval)^2
   ```

   * Evitan divisiones por cero o inestabilidad cuando las varianzas son muy pequeñas.
   * Son los valores clásicos sugeridos en la definición original de SSIM (con coeficientes 0.01 y 0.03).

5. **Fórmula de SSIM**

   ```julia
   num  = (2μx * μy + C1) * (2σxy + C2)
   den  = (μx^2 + μy^2 + C1) * (σx2 + σy2 + C2)

   return num / den
   ```

   Esto implementa:

   $\text{SSIM}(X, Y) =
   \frac{(2\mu_x \mu_y + C_1)(2 \sigma_{xy} + C_2)}
   {(\mu_x^2 + \mu_y^2 + C_1)(\sigma_x^2 + \sigma_y^2 + C_2)}.$

### ¿Por qué es importante?

* A diferencia de MSE/PSNR, SSIM intenta capturar **la percepción visual**:

  * si los patrones estructurales (bordes, texturas) se conservan o no,
  * si el contraste y la luminancia son consistentes.
* Es muy usado en imágenes médicas porque correlaciona mejor con la “calidad visual percibida” que un simple MSE.

En tu proyecto, reportar SSIM junto con PSNR permite mostrar no solo que CS reduce el error numérico, sino que también preserva mejor la **estructura anatómica** que una reconstrucción zero–filled.

---

# Script principal `main_demo.jl` – Pipeline completo de CS en MRI

Este archivo es el **script principal** del proyecto.
Es donde se conectan todos los módulos:

* `Data` → generación/carga de imágenes.
* `Operators` → modelo de adquisición (Fourier + máscara).
* `Masks` → patrón de submuestreo en k-space.
* `ISTA` → algoritmo de reconstrucción CS.
* `Metrics` → evaluación de calidad.

El objetivo de este script es ejecutar el “experimento completo” de *Compressed Sensing* en MRI:

1. Elegir una imagen de referencia $(x_{\text{true}})$ (phantom o imagen real).
2. Simular el k-space completo.
3. Submuestrear con una máscara de densidad variable.
4. Reconstruir:

   * por IFFT zero–filled,
   * por CS (ISTA + wavelets Haar).
5. Calcular métricas (PSNR, SSIM) y tiempos de cómputo.
6. Guardar las imágenes resultantes como PNG.

---

## 1. Importación de paquetes y módulos

```julia
using Images
using FFTW

include("1_data.jl")
include("2_operators.jl")
include("3_masks.jl")
include("4_ista.jl")
include("5_metrics.jl")
```

* `using Images`: para manipular imágenes y guardarlas en PNG (`colorview`, `save`, etc.).

* `using FFTW`: backend de FFT para la transformada de Fourier (aunque la lógica está encapsulada en `Operators.jl`).

* `include(...)` carga los módulos locales:

  * `1_data.jl` → módulo `Data`.
  * `2_operators.jl` → módulo `Operators`.
  * `3_masks.jl` → módulo `Masks`.
  * `4_ista.jl` → módulo `ISTA`.
  * `5_metrics.jl` → módulo `Metrics`.

Esto hace que las definiciones de cada módulo estén disponibles en el espacio de nombres actual.

---

## 2. `using .Data`, `.Operators`, etc.

```julia
using .Data: gaussian_phantom, disk_phantom, shepp_logan_phantom, load_real_image
using .Operators: F_op, F_adj, M_op, M_adj
using .Masks: variable_density_mask
using .ISTA: ista_l1
using .Metrics: psnr, ssim_index
```

* Importa, desde los módulos incluidos:

  * Desde `Data`:

    * `gaussian_phantom`, `disk_phantom`, `shepp_logan_phantom`: phantoms sintéticos.
    * `load_real_image`: carga y normaliza una imagen de MRI real.
  * Desde `Operators`:

    * `F_op`, `F_adj`: transformada de Fourier directa e inversa (modelo MRI).
    * `M_op`, `M_adj`: operador de máscara en k-space.
  * Desde `Masks`:

    * `variable_density_mask`: máscara pseudo–aleatoria para CS.
  * Desde `ISTA`:

    * `ista_l1`: algoritmo de reconstrucción CS con (\ell_1) en wavelets Haar.
  * Desde `Metrics`:

    * `psnr`, `ssim_index`: métricas de calidad.

En resumen: aquí se declara explícitamente qué “bloques” se van a usar para armar el pipeline.

---

## 3. Parámetros globales del experimento

```julia
N = 256
@assert iseven(N) "N debe ser par para la transformada Haar 2D"

accel = 4.0             
center_fraction = 0.10  
λ = 0.02
maxiter = 200
α = 0.8
```

* `N = 256`:

  * Tamaño de la imagen (N×N).
  * Se verifica que sea par porque la transformada Haar 2D requiere dimensiones pares.

* Parámetros de submuestreo:

  * `accel = 4.0`:

    * Aceleración nominal deseada.
    * Aproximadamente se medirá ~1/4 de las líneas periféricas de k-space.
  * `center_fraction = 0.10`:

    * Se asegura que el 10% central de líneas en k-space se mida completo (bajas frecuencias).

* Parámetros de ISTA:

  * `λ = 0.02`: peso de regularización (fuerza de la sparsidad en wavelets).
  * `maxiter = 200`: número máximo de iteraciones.
  * `α = 0.8`: step size del descenso de gradiente.

Estos parámetros controlan qué tan agresivo es el submuestreo y qué tan fuerte es la regularización en la reconstrucción CS.

---

## 4. Elección de la imagen de referencia

```julia
println("Generando imagen de referencia...")
use_real_image = false   # true = usar imagen real, false = usar phantom

if use_real_image
    img_path = "mri_slice.png"
    x_true = load_real_image(img_path; N=N)
else
    # Puedes elegir el phantom que quieras
    # x_true = gaussian_phantom(N)
    # x_true = disk_phantom(N)
    x_true = shepp_logan_phantom(N)
end
```

* Mensaje informativo: inicio de generación de la imagen.
* `use_real_image` permite alternar entre:

  * `true`: usar una imagen real de MRI desde disco.
  * `false`: usar un phantom sintético.

Si `use_real_image == true`:

* `img_path = "mri_slice.png"`: ruta al archivo de imagen.
* `x_true = load_real_image(...)`:

  * Carga la imagen,
  * La pasa a escala de grises,
  * La redimensiona a `N×N`,
  * La normaliza a ([0,1]).

Si `use_real_image == false`:

* Se selecciona un phantom. En este script:

  * Se elige `shepp_logan_phantom(N)` (phantom anatómico tipo cabeza) como referencia principal.

Esta `x_true` es la “verdad terreno” que se utilizará para evaluar la reconstrucción.

---

## 5. Simulación de k-space completo

```julia
println("Calculando k-space completo...")
k_full = F_op(x_true)
```

* Aplica la transformada de Fourier 2D unitaria (centrada) a la imagen verdadera:

$k_{\text{full}} = F x_{\text{true}}.$

* Esto representa la adquisición ideal del resonador **sin submuestreo**:

  * todos los puntos de k-space están disponibles.

---

## 6. Creación de la máscara de submuestreo

```julia
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
```

* `mask = variable_density_mask(...)`:

  * Genera una máscara de densidad variable:

    * bloque central de k-space completo,
    * periferia muestreada de manera pseudo–aleatoria con factor `accel`.

* Se calculan estadísticas del submuestreo:

  * `total_points = N * N`: número total de muestras posibles.
  * `acquired_points = count(!iszero, mask)`: cuántas posiciones en k-space están realmente “encendidas”.
  * `frac_acquired`: fracción de puntos del k-space adquiridos.
  * `R_eff`: factor de aceleración efectivo aproximado.

Se imprimen estos valores para dejar claro **qué tan agresivo** es el submuestreo en el experimento.

---

## 7. Aplicación de la máscara y reconstrucción zero–filled

```julia
println("Aplicando máscara a k-space...")
k_sub = M_op(k_full, mask)

println("Reconstrucción zero-filled (IFFT)...")
t_zf = @elapsed begin
    x_zf = real(F_adj(k_sub))
end
```

* `k_sub = M_op(k_full, mask)`:

  * Aplica la máscara al k-space completo:

    * conserva solo las muestras marcadas por `mask`,
    * las demás se ponen en cero.
  * Simula exactamente la adquisición acelerada.

* Reconstrucción zero–filled:

  * `x_zf = real(F_adj(k_sub))`:

    * aplica IFFT al k-space submuestreado,
    * sin ningún tipo de corrección ni regularización.
  * `@elapsed` mide el tiempo de cómputo `t_zf` para esta operación.

Esta `x_zf` es la referencia clásica “naive” con la que vas a comparar la reconstrucción por CS.

---

## 8. Reconstrucción CS con ISTA + wavelets Haar

```julia
println("Reconstrucción CS con ISTA-L1 (wavelets Haar)...")
t_ista = @elapsed begin
    x_cs = ista_l1(k_sub, mask; λ=λ, maxiter=maxiter, α=α, verbose=true)
end
```

* Llama al algoritmo `ista_l1` pasando:

  * `k_sub`: datos de k-space submuestreado,
  * `mask`: patrón de submuestreo,
  * `λ`, `maxiter`, `α`: parámetros de regularización y optimización.

* Internamente, `ista_l1` resuelve:

  $\min_x ; \tfrac12 | M F x - y |_2^2 + \lambda |W x|_1,$

  con:

  * `F_op`, `F_adj`, `M_op`, `M_adj` (modelo de adquisición),
  * `W_op`, `W_adj` (wavelets Haar 2D),
  * umbralización suave (`soft_threshold`) sobre los coeficientes Haar.

* `t_ista` mide el tiempo total de reconstrucción CS.

El resultado `x_cs` es la imagen reconstruida con *Compressed Sensing*, donde se espera:

* mejor PSNR y SSIM que la zero–filled,
* a costa de un mayor tiempo de cómputo.

---

## 9. Normalización de imágenes para las métricas y visualización

```julia
function norm01(x)
    mn, mx = extrema(x)
    mx == mn && return zeros(size(x))
    return (x .- mn) ./ (mx - mn)
end

x_true_n = norm01(x_true)
x_zf_n   = norm01(x_zf)
x_cs_n   = norm01(x_cs)
```

* Define una función `norm01` que normaliza cualquier imagen a ([0,1]) usando:

  $x_{\text{norm}} = \frac{x - \min(x)}{\max(x) - \min(x)}.$

  Si la imagen es constante (`mx == mn`), devuelve una matriz de ceros para evitar división por cero.

* Aplica `norm01` a:

  * `x_true`,
  * `x_zf`,
  * `x_cs`,
    resultando en `x_true_n`, `x_zf_n`, `x_cs_n`.

**Motivo**:

* Asegurar que las métricas (especialmente PSNR con `maxval=1`) y las imágenes PNG estén en el mismo rango y sean comparables.

---

## 10. Cálculo de métricas PSNR y SSIM

```julia
println("Calculando métricas...")
psnr_zf   = psnr(x_zf_n, x_true_n)
psnr_cs   = psnr(x_cs_n, x_true_n)
ssim_zf_v = ssim_index(x_zf_n, x_true_n)
ssim_cs_v = ssim_index(x_cs_n, x_true_n)
```

* Calcula:

  * `psnr_zf`: PSNR entre zero–filled y la imagen verdadera.
  * `psnr_cs`: PSNR entre CS y la imagen verdadera.
  * `ssim_zf_v`: SSIM entre zero–filled y la verdadera.
  * `ssim_cs_v`: SSIM entre CS y la verdadera.

Esto te da una comparación cuantitativa clara de:

* cuánto mejora CS respecto a zero–filled,
* tanto en error numérico (PSNR) como en similitud estructural (SSIM).

---

## 11. Impresión de resultados y tiempos

```julia
println()
println("=== Resultados ===")
println("Zero-filled:  PSNR = $(round(psnr_zf, digits=2)) dB,  SSIM = $(round(ssim_zf_v, digits=4))")
println("CS (ISTA-L1): PSNR = $(round(psnr_cs, digits=2)) dB,  SSIM = $(round(ssim_cs_v, digits=4))")

println()
println("=== Tiempos de cómputo ===")
println("Zero-filled (IFFT): $(round(t_zf, digits=4)) s")
println("CS (ISTA + wavelets): $(round(t_ista, digits=4)) s")
```

* Imprime en consola una tabla textual con:

  * PSNR y SSIM para zero–filled,
  * PSNR y SSIM para CS.

* Luego, tiempos:

  * tiempo de reconstrucción zero–filled,
  * tiempo de reconstrucción CS.

Esto resume el trade–off fundamental del proyecto:

* **Menos datos** + CS
  vs
* Reconstrucción simple zero–filled,

mostrando explícitamente que CS:

* **mejora la calidad de imagen**,
* pero **requiere más tiempo de cómputo**.

---

## 12. Guardado de imágenes PNG

```julia
println("Guardando imágenes PNG...")
save("x_true.png", colorview(Gray, x_true_n))
save("x_zf.png",   colorview(Gray, x_zf_n))
save("x_cs.png",   colorview(Gray, x_cs_n))

println("Listo. Imágenes guardadas en el directorio actual.")
```

* Guarda las tres imágenes normalizadas como:

  * `x_true.png` → imagen de referencia.
  * `x_zf.png` → reconstrucción zero–filled.
  * `x_cs.png` → reconstrucción CS.

* `colorview(Gray, x_*)` las convierte a imagen en escala de grises compatible con `save`.

Estas imágenes son las que después puedes poner en el informe para:

* mostrar visualmente el efecto del submuestreo,
* comparar artefactos de aliasing vs. resultado CS,
* ilustrar las diferencias que justifican los números de PSNR y SSIM.

---

## 13. Resumen conceptual del script

Este script implementa **todo el pipeline práctico** de tu solución:

1. **Definir el problema**:

   * Imagen de referencia $(x_{\text{true}})$,
   * modelo de adquisición $(y = MFx)$,
   * máscara de submuestreo de densidad variable (CS).

2. **Simular la adquisición acelerada**:

   * calcular k-space completo,
   * aplicar la máscara,
   * obtener los datos reducidos `k_sub`.

3. **Comparar métodos de reconstrucción**:

   * Zero–filled: IFFT directa,
   * CS con ISTA + wavelets Haar.

4. **Medir y documentar el desempeño**:

   * PSNR y SSIM para cada método,
   * tiempos de cómputo,
   * guardar imágenes para análisis visual.