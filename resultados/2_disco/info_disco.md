# Resultado – CS con *phantom* disco

Este experimento corresponde a la ejecución de `6_main_demo.jl` usando como imagen de referencia el **phantom tipo disco** (`disk_phantom(N)`).

Las imágenes asociadas son:

* `x_true.png` – disco ideal, homogéneo y bien definido.
* `x_zf.png` – reconstrucción *zero–filled* (IFFT de k-space submuestreado).
* `x_cs.png` – reconstrucción por *Compressed Sensing* (ISTA + wavelets Haar).

---

## Parámetros de adquisición simulada

* Tamaño de imagen: $N = 256 \Rightarrow 256 \times 256 = 65{,}536$ puntos de k-space.
* Máscara de submuestreo de densidad variable con:

  * `accel = 4.0`
  * `center_fraction = 0.10`

Con esta máscara, el script reporta:

* Puntos totales en k-space: **65.536**
* Puntos adquiridos: **21.504**
* Fracción adquirida: **32,81 %**
* Factor de aceleración efectivo: $R_{\text{eff}} \approx 3{,}05$
  
Es decir, se simula reducir el tiempo de adquisición a ~1/3 respecto a muestrear todo k-space.

---

## Reconstrucciones y comparación visual

* **Referencia (`x_true`)**:
  Disco blanco uniforme sobre fondo negro, bordes nítidos.

* **Zero–filled (`x_zf`)**:

  * Bordes del disco con artefactos de oscilación (anillos verticales y bandas alrededor).
  * Pérdida de homogeneidad dentro del disco.
  * Representa lo que pasa si solo se hace IFFT con k-space submuestreado.

* **CS (ISTA + Haar, `x_cs`)**:

  * Bordes mucho más limpios que en la reconstrucción zero–filled.
  * Interior del disco más homogéneo; se reducen notablemente los artefactos.
  * Queda todavía una textura muy sutil, pero la forma del disco se aproxima muy bien a la referencia.

---

## Métricas y tiempos

Resultados que entrega el script:

```text
=== Resultados ===
Zero-filled:  PSNR = 21.17 dB,  SSIM = 0.6734
CS (ISTA-L1): PSNR = 43.59 dB,  SSIM = 0.9969

=== Tiempos de cómputo ===
Zero-filled (IFFT): 0.006 s
CS (ISTA + wavelets): 2.7486 s
```

Interpretación rápida:

* **CS mejora mucho la calidad** respecto a zero–filled:

  * PSNR sube de ~21 dB a ~43,6 dB.
  * SSIM pasa de ~0,67 a ~0,997 (casi idéntico a la referencia).
* **Costo computacional**:

  * IFFT zero–filled es casi instantánea (~0,006 s).
  * Reconstrucción CS tarda unos ~2,75 s debido al algoritmo iterativo.

En resumen, con solo ~33 % de los puntos de k-space el método CS recupera un disco prácticamente perfecto, mientras que la reconstrucción zero–filled muestra artefactos claros en los bordes.
