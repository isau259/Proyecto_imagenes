# Resultado – CS con *phantom* gaussiano

Este experimento corresponde a la ejecución de `6_main_demo.jl` usando como imagen de referencia un **phantom gaussiano 2D**.

## Parámetros principales

* Tamaño de imagen: $N = 256 \Rightarrow 256 \times 256 = 65{,}536$ píxeles / puntos de k-space.
* Imagen de referencia: `gaussian_phantom(N)` → “spot” gaussiano suave en el centro del FOV.
* Submuestreo en k-space:

  * `accel = 4.0`
  * `center_fraction = 0.10`

A partir de la máscara generada:

* Puntos totales en k-space: **65.536**
* Puntos adquiridos: **21.504**
* Fracción adquirida: **32,81 %**
* Factor de aceleración efectivo: $R_{\text{eff}} \approx 3{,}05$

Interpretación: se emula una reducción del tiempo de adquisición de ~3× manteniendo solo un tercio de los puntos de k-space.

## Reconstrucciones

Se obtienen tres imágenes:

* `x_true.png` – phantom gaussiano original.
* `x_zf.png` – reconstrucción *zero-filled* (IFFT directa de k-space submuestreado).
* `x_cs.png` – reconstrucción por *Compressed Sensing* (ISTA + wavelets Haar 2D).

Visualmente:

* Las tres imágenes (`x_true`, `x_zf`, `x_cs`) son prácticamente indistinguibles: un punto gaussiano suave sobre fondo negro.
* No se aprecian artefactos visibles ni pérdida de contraste en este caso tan simple.

## Métricas y tiempos

Resultados del script:

```text
=== Resultados ===
Zero-filled:  PSNR = 182.24 dB,  SSIM = 1.0
CS (ISTA-L1): PSNR = 46.7 dB,  SSIM = 0.9984

=== Tiempos de cómputo ===
Zero-filled (IFFT): 0.0136 s
CS (ISTA + wavelets): 13.567 s
```

Comentarios:

* La métrica SSIM es ≈1 para ambas reconstrucciones, coherente con el hecho de que el phantom gaussiano es muy simple y el submuestreo utilizado no destruye información relevante.
* El PSNR de CS es alto (≈46,7 dB) y el de zero-filled aparece saturado numéricamente en esta corrida; en cualquier caso, **no hay diferencia visual apreciable** entre las tres imágenes.
* El coste computacional sí cambia: la IFFT *zero-filled* es casi instantánea (~0,014 s), mientras que la reconstrucción CS requiere varios segundos (~13,6 s) debido al proceso iterativo con wavelets.

En resumen, este experimento con phantom gaussiano sirve como caso de prueba básico: k-space se submuestrea hasta un ~33 %, pero, dada la simplicidad espacial de la imagen, tanto la reconstrucción zero-filled como la CS recuperan el phantom sin artefactos visibles.
