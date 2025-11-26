# Resultado – CS con imagen real de cerebro (MRI T2)

Este experimento corresponde a la ejecución de `6_main_demo.jl` usando como imagen de referencia una **corte axial real de cerebro** (MRI ponderada en T2), cargada desde archivo y reescalada a (256 \times 256).

Las imágenes asociadas son:

* `x_true.png` – imagen MRI original (normalizada, sin ruido añadido).
* `x_zf.png` – reconstrucción *zero–filled* (IFFT directa de k-space submuestreado).
* `x_cs.png` – reconstrucción por *Compressed Sensing* (ISTA + wavelets Haar).

---

## Parámetros de adquisición simulada

* Tamaño de imagen: $N = 256 \Rightarrow 65{,}536$ puntos de k-space.
* Máscara de submuestreo de densidad variable:

  * `accel = 4.0`
  * `center_fraction = 0.10`

Resultados de la máscara:

* Puntos totales en k-space: **65.536**
* Puntos adquiridos: **21.504**
* Fracción adquirida: **32,81 %**
* Factor de aceleración efectivo: $R_{\text{eff}} \approx 3{,}05$

Interpretación: se emula adquirir solo ~1/3 de las líneas de k-space, equivalente a reducir el tiempo de adquisición en un factor ~3.

---

## Comparación visual

* **Referencia (`x_true`)**
  Imagen cerebral nítida, con buena diferenciación de sustancia gris/blanca y LCR, y contornos bien definidos.

* **Zero–filled (`x_zf`)**

  * Se observan **aliasing y patrones periódicos** que atraviesan el parénquima.
  * Hay pérdida de contraste local y textura “lavada” en varias regiones.
  * Aun así, la anatomía general sigue reconocible, pero la calidad se ve claramente degradada.

* **CS (ISTA + Haar, `x_cs`)**

  * Reduce buena parte de los artefactos oscilatorios visibles en la reconstrucción zero–filled.
  * Mejora el contraste de bordes y la definición de surcos y ventrículos.
  * Persisten artefactos residuales y cierta suavización de detalles finos, lo que refleja que el modelo de sparsidad en wavelets es simple frente a la complejidad real de la anatomía.

---

## Métricas y tiempos

Salida del script:

```text
=== Resultados ===
Zero-filled:  PSNR = 17.24 dB,  SSIM = 0.7422
CS (ISTA-L1): PSNR = 21.62 dB,  SSIM = 0.8674

=== Tiempos de cómputo ===
Zero-filled (IFFT): 0.0097 s
CS (ISTA + wavelets): 6.445 s
```

Interpretación:

* **Calidad de imagen**

  * CS mejora el PSNR en ~4,4 dB y el SSIM de ~0,74 a ~0,87.
  * La mejora es clara, pero no tan extrema como en los phantoms simples: la imagen real tiene estructuras mucho más ricas, por lo que el modelo de sparsidad no logra recuperar todos los detalles.

* **Costo computacional**

  * Reconstrucción zero–filled: ~0,01 s (IFFT directa).
  * Reconstrucción CS: ~6,45 s (algoritmo iterativo ISTA + wavelets).

En resumen, en una **imagen MRI real de cerebro**, el submuestreo al 33 % genera aliasing significativo si solo se aplica IFFT; al aplicar **Compressed Sensing**, se logra una imagen visualmente más cercana a la original y con métricas objetivas mejores, a costa de un tiempo de cómputo varias órdenes de magnitud mayor. Este caso ilustra el compromiso real entre calidad, aceleración y complejidad del modelo de reconstrucción.
