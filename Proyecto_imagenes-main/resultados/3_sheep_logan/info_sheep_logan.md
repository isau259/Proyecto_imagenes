# Resultado – CS con *phantom* Shepp–Logan

Este experimento corresponde a la ejecución de `6_main_demo.jl` usando como imagen de referencia el **phantom de Shepp–Logan** (`shepp_logan_phantom(N)`), que modela una sección axial de cabeza con varias estructuras internas (lóbulos, masa central, etc.).

Las imágenes asociadas son:

* `x_true.png` – phantom Shepp–Logan original.
* `x_zf.png` – reconstrucción *zero–filled* (IFFT de k-space submuestreado).
* `x_cs.png` – reconstrucción por *Compressed Sensing* (ISTA + wavelets Haar).

---

## Parámetros de adquisición simulada

* Tamaño de imagen: $N = 256 \Rightarrow 256 \times 256 = 65{,}536$ puntos de k-space.
* Máscara de submuestreo de densidad variable:

  * `accel = 4.0`
  * `center_fraction = 0.10`

El script reporta:

* Puntos totales en k-space: **65.536**
* Puntos adquiridos: **21.504**
* Fracción adquirida: **32,81 %**
* Factor de aceleración efectivo: $R_{\text{eff}} \approx 3{,}05$

Se emula, por tanto, una reducción del tiempo de adquisición a ~1/3 respecto al muestreo completo.

---

## Reconstrucciones y comparación visual

* **Referencia (`x_true`)**
  Phantom Shepp–Logan con contorno elíptico blanco y estructuras internas bien definidas en diferentes niveles de gris.

* **Zero–filled (`x_zf`)**

  * Fuerte presencia de artefactos de aliasing, en forma de anillos y patrones oscilatorios horizontales.
  * Las estructuras internas se ven “lavadas” y con bajo contraste.
  * Es un ejemplo claro de cómo el submuestreo directo degrada la imagen en un objeto complejo.

* **CS (ISTA + Haar, `x_cs`)**

  * Recupera mucho mejor la forma global del cráneo y las cavidades internas.
  * Se reducen notablemente los patrones de aliasing presentes en `x_zf`.
  * Todavía se observan artefactos residuales y cierta pérdida de nitidez en detalles finos, pero la imagen es claramente más diagnóstica que la zero–filled.

---

## Métricas y tiempos

Del output del script:

```text
=== Resultados ===
Zero-filled:  PSNR = 16.62 dB,  SSIM = 0.6932
CS (ISTA-L1): PSNR = 23.85 dB,  SSIM = 0.9237

=== Tiempos de cómputo ===
Zero-filled (IFFT): 0.4038 s
CS (ISTA + wavelets): 5.4365 s
```

Interpretación rápida:

* **Calidad de imagen**

  * CS mejora el PSNR en ~7 dB (de 16,6 a 23,9 dB).
  * SSIM sube de ~0,69 (calidad moderada, con aliasing visible) a ~0,92 (estructura mucho más parecida al phantom original).

* **Costo computacional**

  * La IFFT zero–filled toma ~0,40 s.
  * La reconstrucción CS con ISTA + wavelets toma ~5,44 s.

En resumen, para este phantom más complejo, el submuestreo del 67 % de k-space introduce aliasing severo en la reconstrucción zero–filled, mientras que **Compressed Sensing** logra recuperar de forma mucho más fiel la anatomía simulada, a costa de un mayor tiempo de cómputo.
