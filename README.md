# 💻 Simulasi Kinetik Jaringan Metabolik dengan Allosteric Feedback Inhibition

## 📌 Gambaran Sistem

Sistem metabolik yang disimulasikan terdiri atas satu jalur utama pembentukan produk dan satu jalur samping (byproduct). Produk akhir (**P**) memberikan umpan balik negatif (*allosteric feedback inhibition*) terhadap reaksi pertama (**v₁**) sehingga laju pembentukan metabolit awal akan menurun seiring meningkatnya konsentrasi produk.

```text
                    Allosteric Inhibition (non-competitive)
         ┌──────────────────────────────────────────────┐
         │                                              ▼
  X ──[v1]──► A ──[v2]──► B ──[v3]──► P
                   │
                  [v4]
                   │
                   ▼
               Byproduct
```

---

## ⚙️ Parameter Simulasi

| Parameter | Keterangan                     | Nilai |
| --------- | ------------------------------ | ----- |
| V₁,max    | Laju maksimum reaksi v₁        | 5.0   |
| Kₘ₁       | Michaelis constant             | 2.0   |
| Kᵢ        | Konstanta inhibisi             | 3.0   |
| X         | Konsentrasi substrat eksternal | 10    |
| k₂        | Konstanta laju A → B           | 1.0   |
| k₃        | Konstanta laju B → P           | 0.8   |
| k₄        | Konstanta laju A → Byproduct   | 0.3   |

---

## 🧮 Persamaan Kinetik

### Reaksi v₁ (Michaelis-Menten dengan Inhibisi Non-Kompetitif)

[
v_1=\frac{V_{1,max}[X]}
{(K_{m1}+[X])\left(1+\frac{[P]}{K_i}\right)}
]

### Reaksi v₂, v₃, dan v₄

[
v_2=k_2[A]
]

[
v_3=k_3[B]
]

[
v_4=k_4[A]
]

---

## 📈 Sistem Persamaan Diferensial

### Metabolit A

[
\frac{d[A]}{dt}=v_1-v_2-v_4
]

### Metabolit B

[
\frac{d[B]}{dt}=v_2-v_3
]

### Produk P

[
\frac{d[P]}{dt}=v_3
]

---

## 🐍 Script Python Simulasi

```python
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parameter
V1max = 5.0
Km1 = 2.0
Ki = 3.0
X = 10.0
k2 = 1.0
k3 = 0.8
k4 = 0.3

# Sistem ODE
def metabolic_system(t, y):
    A, B, P = y

    v1 = (V1max * X) / ((Km1 + X) * (1 + P / Ki))
    v2 = k2 * A
    v3 = k3 * B
    v4 = k4 * A

    dA_dt = v1 - v2 - v4
    dB_dt = v2 - v3
    dP_dt = v3

    return [dA_dt, dB_dt, dP_dt]

# Kondisi awal
y0 = [0, 0, 0]

# Waktu simulasi
t_span = (0, 50)
t_eval = np.linspace(0, 50, 1000)

# Penyelesaian numerik
solution = solve_ivp(
    metabolic_system,
    t_span,
    y0,
    t_eval=t_eval,
    method='RK45'
)

# Visualisasi
plt.figure(figsize=(8,5))

plt.plot(solution.t, solution.y[0], label='[A]')
plt.plot(solution.t, solution.y[1], label='[B]')
plt.plot(solution.t, solution.y[2], label='[P]')

plt.xlabel('Waktu')
plt.ylabel('Konsentrasi')
plt.title('Simulasi Kinetik Metabolik')
plt.legend()
plt.grid(True)

plt.show()
```

---

## 📊 Hasil Simulasi

Simpan hasil visualisasi sebagai:

```text
simulasi_metabolik.png
```

Lalu tampilkan pada README menggunakan:

```markdown
![Simulasi Kinetik](simulasi_metabolik.png)
```

---

## 🔍 Interpretasi Hasil

* Konsentrasi **A** meningkat pada fase awal karena reaksi v₁ masih bekerja maksimal.
* Konsentrasi **B** meningkat setelah terbentuknya A dan kemudian mendekati kondisi tunak (*steady state*).
* Produk **P** terus terakumulasi sepanjang simulasi.
* Akumulasi **P** menghambat reaksi v₁ melalui mekanisme *allosteric feedback inhibition*.
* Semakin tinggi konsentrasi **P**, semakin rendah laju pembentukan A sehingga keseluruhan flux metabolik menurun.
* Mekanisme ini merepresentasikan sistem regulasi biologis yang menjaga kestabilan metabolisme sel.

---

## 📚 Referensi

1. Palsson, B.Ø. (2015). *Systems Biology: Constraint-based Reconstruction and Analysis*. Cambridge University Press.
2. Fell, D.A. (1997). *Understanding the Control of Metabolism*. Portland Press.
3. Strogatz, S.H. (2018). *Nonlinear Dynamics and Chaos*. CRC Press.
4. Virtanen, P. et al. (2020). *SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python*. Nature Methods, 17, 261–272.
