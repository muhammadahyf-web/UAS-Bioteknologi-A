# 📘 UAS Bioteknologi — Systems Biology
### Matakuliah: BISB211605 | Fakultas Biologi, Universitas Gadjah Mada
### Dosen: Matin Nuhamunada, M.Sc., Ph.D. | Deadline: 19 Juni 2026 Pukul 23.45

---

## 📋 Daftar Isi

1. [Pemahaman Skenario dan Sistem](#-skenario-sistem-metabolik)
2. [Question 1 — Memilih Modeling Framework](#-question-1--memilih-modeling-framework)
3. [Question 2A — Stoichiometric Matrix](#-question-2a--stoichiometric-matrix)
4. [Question 2B — Persamaan Diferensial (ODE)](#-question-2b--kinetic-differential-equations)
5. [Question 2C — Simulasi Kinetik (Bonus)](#-question-2c-bonus--simulasi-kinetik-dengan-allosteric-inhibition)
6. [Referensi](#-referensi)

---

## 🧬 Skenario Sistem Metabolik

### Gambaran Umum

Kita berperan sebagai **metabolic engineer** yang mengoptimalkan mikroorganisme untuk memproduksi suatu biochemical bernilai tinggi (disebut **Product P**). Di dalam sel, terdapat jaringan reaksi metabolik berikut:

```
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

### Komponen Sistem

| Simbol | Jenis | Keterangan |
|--------|-------|------------|
| **X** | Substrat eksternal | Bahan baku awal; dijaga konstan dari luar sel |
| **A** | Metabolit internal | Produk antara pertama; titik percabangan (branch point) |
| **B** | Metabolit internal | Produk antara kedua |
| **P** | Metabolit internal | Produk akhir yang diinginkan |
| **Byproduct** | Keluaran sampingan | Produk tidak diinginkan yang membuang flux dari jalur utama |

### Reaksi dalam Sistem

| Reaksi | Transformasi | Jenis Kinetik | Keterangan |
|--------|-------------|---------------|------------|
| **v₁** | X → A | Michaelis-Menten | Langkah komitmen awal; dihambat secara alosterik oleh P |
| **v₂** | A → B | First-order (k₂) | Konversi intermediat |
| **v₃** | B → P | First-order (k₃) | Langkah produksi akhir |
| **v₄** | A → Byproduct | First-order (k₄) | Kebocoran flux ke jalur sampingan |

### Karakteristik Kunci Sistem

- Substrat **X** bersifat eksternal dan **konsentrasinya dijaga konstan** (tidak termasuk dalam ODE)
- **Allosteric feedback inhibition**: P menghambat v₁ secara non-kompetitif — semakin banyak P terbentuk, semakin lambat v₁ bekerja
- **Branch point** terjadi di metabolit A: flux terbagi antara jalur produktif (v₂ → B → P) dan jalur buangan (v₄ → Byproduct)

---

## 🔬 Question 1 — Memilih Modeling Framework

### Latar Belakang: Tiga Framework yang Tersedia

Sebelum menjawab, penting untuk memahami kapabilitas dan batasan masing-masing framework:

#### A. Kinetic Modeling (ODE-based)
Framework ini memodelkan sistem metabolik sebagai sekumpulan **persamaan diferensial biasa (Ordinary Differential Equations)**. Untuk setiap metabolit internal, ditulis persamaan yang menggambarkan laju perubahan konsentrasinya terhadap waktu berdasarkan hukum kekekalan massa. Keunggulan utama ODE adalah kemampuannya menangkap **dinamika temporal** — yaitu bagaimana konsentrasi metabolit berubah dari detik ke detik. Namun, ODE membutuhkan **parameter kinetik lengkap** untuk setiap reaksi, yang seringkali sulit diperoleh untuk jaringan skala besar.

#### B. Constraint-Based Modeling (Flux Balance Analysis / FBA)
FBA menggunakan asumsi **steady-state**, bahwa laju perubahan konsentrasi semua metabolit internal adalah nol pada kondisi operasional, sehingga berlaku:

```
S · v = 0
```

di mana **S** adalah stoichiometric matrix dan **v** adalah vektor flux reaksi. Dengan mendefinisikan fungsi objektif dan constraint kapasitas flux, FBA menggunakan **linear programming** untuk menemukan distribusi flux optimal. Keunggulan besarnya adalah kemampuan menangani jaringan metabolik **skala genom** (ribuan reaksi) tanpa memerlukan data kinetik.

#### C. Network Topology / Graph Theory
Pendekatan ini memodelkan jaringan metabolisme sebagai **graf berarah (directed graph)** di mana metabolit adalah node dan reaksi adalah edge. Analisis topologi seperti degree centrality dan betweenness centrality digunakan untuk memahami **struktur arsitektural** jaringan — bukan untuk menghitung flux atau dinamika konsentrasi.

---

### Jawaban Q1

#### Skenario 1 → **B: Constraint-Based Modeling (FBA)**

> *"Menghitung maximum theoretical yield produk di seluruh jaringan metabolik skala genom"*

FBA adalah pilihan yang paling tepat karena pertanyaan yang diajukan adalah kuantitatif dalam konteks **optimasi di jaringan skala besar**. FBA bekerja dengan menyelesaikan permasalahan optimasi linear — memaksimalkan fungsi objektif (produksi P) dengan constraint S·v = 0 dan batas kapasitas flux. Fokus FBA tepat sasaran karena menghasilkan **yield maksimum teoritis** dan **produksi maksimum** pada **genome-scale network** tanpa memerlukan data kinetik untuk setiap reaksi. Kinetic Modeling (ODE) tidak relevan di sini karena membutuhkan parameter kinetik yang tidak mungkin tersedia untuk ribuan reaksi sekaligus, dan Graph Theory hanya bisa menjawab pertanyaan struktural, bukan kuantitatif tentang yield.

#### Skenario 2 → **A: Kinetic Modeling (ODE-based)**

> *"Mensimulasi akumulasi metabolit B yang toksik secara real-time selama 48 jam fermentasi"*

Kinetic Modeling berbasis ODE adalah satu-satunya framework yang mampu menjawab pertanyaan ini karena mengandung **dimensi waktu yang eksplisit**. ODE secara langsung memodelkan laju perubahan konsentrasi setiap metabolit sebagai fungsi waktu, sehingga dapat memprediksi **kapan tepatnya** [B] mencapai ambang toksik. Informasi ini krusial untuk merancang strategi intervensi fermentasi yang tepat waktu. FBA tidak dapat digunakan karena ia mengasumsikan steady-state dan tidak mengenal dimensi waktu sama sekali, sedangkan Graph Theory hanya menganalisis struktur jaringan tanpa informasi kuantitatif kinetik.

#### Skenario 3 → **C: Network Topology / Graph Theory**

> *"Mengidentifikasi structural hub di organisme baru agar gen knockout tidak mengganggu jalur esensial"*

Graph Theory adalah pilihan yang paling tepat karena pertanyaan ini murni bersifat **struktural** — kita ingin mengetahui arsitektur konektivitas jaringan, bukan dinamika atau optimasi flux. Dengan merepresentasikan metabolisme sebagai graf berarah, kita dapat menghitung metrik topologi seperti **degree centrality** (berapa banyak reaksi yang terhubung ke suatu metabolit) dan **betweenness centrality** (seberapa sering metabolit tersebut menjadi jembatan antar bagian jaringan). Metabolit dengan nilai centrality tinggi adalah hub yang jika terganggu dapat menyebabkan kegagalan sistemik. Penting dicatat bahwa pada organisme yang baru disequens, data kinetik belum tentu tersedia, sehingga ODE tidak praktis digunakan.

---

## 📊 Question 2A — Stoichiometric Matrix

### Konsep Dasar

Stoichiometric Matrix (S) adalah representasi matematis kompak dari seluruh reaksi dalam jaringan metabolik. Matriks ini berfungsi sebagai fondasi matematis untuk analisis FBA maupun penulisan ODE. Aturan pengisian:

| Nilai | Makna |
|-------|-------|
| **+1** | Metabolit **diproduksi** oleh reaksi tersebut |
| **-1** | Metabolit **dikonsumsi** oleh reaksi tersebut |
| **0** | Metabolit **tidak terlibat** dalam reaksi tersebut |

> **Catatan:** Hanya metabolit **internal** (A, B, P) yang menjadi baris matriks. X adalah substrat eksternal (konstan), dan Byproduct dianggap keluar dari sistem, sehingga keduanya tidak dimasukkan.

### Analisis Setiap Reaksi

| Reaksi | Metabolit A | Metabolit B | Product P | Alasan |
|--------|:-----------:|:-----------:|:---------:|--------|
| **v₁** (X→A) | +1 | 0 | 0 | A diproduksi; B dan P tidak terlibat |
| **v₂** (A→B) | -1 | +1 | 0 | A dikonsumsi, B diproduksi; P tidak terlibat |
| **v₃** (B→P) | 0 | -1 | +1 | B dikonsumsi, P diproduksi; A tidak terlibat |
| **v₄** (A→Byp) | -1 | 0 | 0 | A dikonsumsi; B dan P tidak terlibat |

### Jawaban Q2A — Stoichiometric Matrix Final

```
        v₁   v₂   v₃   v₄
A   [  +1   -1    0   -1  ]
B   [   0   +1   -1    0  ]
P   [   0    0   +1    0  ]
```

Matriks ini berukuran **3×4** (3 metabolit internal × 4 reaksi). Setiap baris matriks secara langsung merepresentasikan persamaan keseimbangan massa untuk metabolit yang bersangkutan. Pada kondisi steady-state dalam analisis FBA, persamaan S·v = 0 harus terpenuhi untuk setiap baris.

---

## ⚗️ Question 2B — Kinetic Differential Equations

### Konsep Dasar

ODE untuk sistem metabolik diturunkan dari **hukum kekekalan massa**:

```
d[metabolit]/dt = Σ flux_produksi − Σ flux_konsumsi
```

Persamaan ini identik dengan baris yang bersesuaian pada stoichiometric matrix dikalikan dengan vektor flux.

### Ekspresi Kinetik yang Digunakan

**Reaksi v₁ — Michaelis-Menten Kinetics** *(abaikan inhibisi allosterik untuk Q2B sesuai hint soal)*

Reaksi v₁ dikatalisis oleh enzim yang memiliki kapasitas maksimal (Vmax) dan afinitas substrat tertentu (Km). Ketika konsentrasi substrat jauh melebihi Km, laju mendekati Vmax (saturasi):

```
v₁ = (Vmax · [X]) / (Km1 + [X])
```

**Reaksi v₂, v₃, v₄ — First-Order Kinetics**

Laju reaksi berbanding lurus dengan konsentrasi substratnya:

```
v₂ = k₂ · [A]
v₃ = k₃ · [B]
v₄ = k₄ · [A]
```

### Jawaban Q2B — Sistem Persamaan ODE

**Persamaan 1: Metabolit A**

A diproduksi oleh v₁ (dari X) dan dikonsumsi oleh v₂ (menuju B) serta v₄ (menuju Byproduct):

```
d[A]/dt = v₁ - v₂ - v₄ = (Vmax·[X])/(Km1+[X]) - k₂[A] - k₄[A]
```

**Persamaan 2: Metabolit B**

B hanya diproduksi oleh v₂ (dari A) dan hanya dikonsumsi oleh v₃ (menuju P):

```
d[B]/dt = v₂ - v₃ = k₂[A] - k₃[B]
```

**Persamaan 3: Product P**

P hanya diproduksi oleh v₃ (dari B) dan tidak dikonsumsi oleh reaksi manapun dalam sistem ini:

```
d[P]/dt = v₃ = k₃[B]
```

### Interpretasi Sistem ODE

Ketiga persamaan ini membentuk sistem yang **saling bergantung (coupled)**: konsentrasi A mempengaruhi laju pembentukan B (melalui v₂) dan pembentukan Byproduct (melalui v₄), dan konsentrasi B mempengaruhi laju pembentukan P (melalui v₃). Oleh karena itu, sistem ini tidak dapat diselesaikan secara analitik satu per satu, melainkan harus diselesaikan secara **simultan dan numerik** menggunakan metode seperti Runge-Kutta (RK45).

---

## 💻 Question 2C (Bonus) — Simulasi Kinetik dengan Allosteric Inhibition

### Penjelasan Inhibisi Non-Kompetitif

Pada bagian bonus ini, model dilengkapi dengan **efek inhibisi alosterik non-kompetitif** dari Product P terhadap enzim yang mengkatalisis v₁. Dalam inhibisi non-kompetitif, inhibitor (P) tidak bersaing dengan substrat (X) pada sisi aktif enzim, melainkan berikatan pada **sisi alosterik** yang terpisah. Akibatnya, Km tidak berubah (afinitas enzim terhadap substrat tetap sama), namun Vmax efektif menurun (kapasitas maksimal enzim berkurang). Secara matematis, inhibisi ini direpresentasikan dengan menambahkan **faktor inhibisi** pada penyebut persamaan Michaelis-Menten:

```
v₁ = (V1max · [X]) / ((Km1 + [X]) · (1 + [P]/Ki))
```

Ketika [P] = 0, faktor inhibisi bernilai 1 sehingga persamaan kembali ke bentuk Michaelis-Menten standar. Seiring meningkatnya [P], penyebut membesar sehingga v₁ semakin kecil. Mekanisme ini adalah **regulasi homeostasis alami** sel yang mencegah overakumulasi produk.

### Parameter Simulasi

| Parameter | Makna | Nilai |
|-----------|-------|-------|
| V₁,max | Laju maksimal reaksi v₁ | 5.0 |
| Km₁ | Michaelis constant untuk v₁ | 2.0 |
| Kᵢ | Konstanta inhibisi | 3.0 |
| X | Konsentrasi substrat eksternal (konstan) | 10 |
| k₂ | Rate constant reaksi v₂ (A → B) | 1.0 |
| k₃ | Rate constant reaksi v₃ (B → P) | 0.8 |
| k₄ | Rate constant reaksi v₄ (A → Byproduct) | 0.3 |

### Script Python Simulasi

```python
# ============================================================
# SIMULASI KINETIK JARINGAN METABOLIK
# Matakuliah Bioteknologi (BISB211605)
# Fakultas Biologi, Universitas Gadjah Mada
# ============================================================
# Referensi: https://github.com/lab-biotek-bio-ugm/S1_BISB211605_Biotechnology
# ============================================================

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ─────────────────────────────────────────────
# 1. DEFINISI PARAMETER
# ─────────────────────────────────────────────
V1max = 5.0    # Laju maksimal reaksi v1 (unit/waktu)
Km1   = 2.0    # Michaelis constant untuk v1
Ki    = 3.0    # Konstanta inhibisi non-kompetitif
X     = 10.0   # Konsentrasi substrat eksternal (konstan)
k2    = 1.0    # Rate constant v2: A → B  (1/waktu)
k3    = 0.8    # Rate constant v3: B → P  (1/waktu)
k4    = 0.3    # Rate constant v4: A → Byproduct (1/waktu)


# ─────────────────────────────────────────────
# 2. DEFINISI SISTEM ODE
# ─────────────────────────────────────────────
def sistem_ode_dengan_inhibisi(t, y):
    """
    Sistem ODE dengan allosteric feedback inhibition.
    y = [A, B, P]
    """
    A, B, P = y

    # v1: Michaelis-Menten dengan non-competitive inhibition oleh P
    v1 = (V1max * X) / ((Km1 + X) * (1 + P / Ki))

    # v2, v3, v4: First-order kinetics
    v2 = k2 * A
    v3 = k3 * B
    v4 = k4 * A

    # Persamaan ODE (hukum kekekalan massa)
    dA_dt = v1 - v2 - v4   # A: diproduksi v1, dikonsumsi v2 dan v4
    dB_dt = v2 - v3         # B: diproduksi v2, dikonsumsi v3
    dP_dt = v3              # P: hanya diproduksi v3 (tidak dikonsumsi)

    return [dA_dt, dB_dt, dP_dt]


def sistem_ode_tanpa_inhibisi(t, y):
    """
    Sistem ODE TANPA allosteric feedback inhibition (baseline).
    Digunakan untuk perbandingan.
    """
    A, B, P = y

    # v1: Michaelis-Menten standar (tanpa inhibisi)
    v1 = (V1max * X) / (Km1 + X)

    v2 = k2 * A
    v3 = k3 * B
    v4 = k4 * A

    dA_dt = v1 - v2 - v4
    dB_dt = v2 - v3
    dP_dt = v3

    return [dA_dt, dB_dt, dP_dt]


# ─────────────────────────────────────────────
# 3. KONDISI AWAL DAN WAKTU SIMULASI
# ─────────────────────────────────────────────
y0     = [0.0, 0.0, 0.0]            # A=0, B=0, P=0 pada t=0
t_span = (0, 50)                     # Simulasi dari t=0 sampai t=50
t_eval = np.linspace(0, 50, 1000)   # 1000 titik waktu


# ─────────────────────────────────────────────
# 4. SELESAIKAN SISTEM ODE SECARA NUMERIK
# ─────────────────────────────────────────────
sol_inhibisi = solve_ivp(sistem_ode_dengan_inhibisi, t_span, y0,
                         t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)

sol_baseline = solve_ivp(sistem_ode_tanpa_inhibisi, t_span, y0,
                         t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)


# ─────────────────────────────────────────────
# 5. HITUNG FLUX v1 SEPANJANG WAKTU
# ─────────────────────────────────────────────
P_inh   = sol_inhibisi.y[2]
v1_inh  = (V1max * X) / ((Km1 + X) * (1 + P_inh / Ki))
v1_base = np.full_like(v1_inh, (V1max * X) / (Km1 + X))


# ─────────────────────────────────────────────
# 6. VISUALISASI HASIL SIMULASI
# ─────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

warna = {
    'A'   : '#2196F3',
    'B'   : '#FF9800',
    'P'   : '#4CAF50',
    'v1'  : '#E91E63',
    'base': '#9E9E9E'
}

# Panel 1: Dinamika konsentrasi dengan inhibisi
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(sol_inhibisi.t, sol_inhibisi.y[0], color=warna['A'], lw=2.5, label='[A] Intermediat A')
ax1.plot(sol_inhibisi.t, sol_inhibisi.y[1], color=warna['B'], lw=2.5, label='[B] Intermediat B')
ax1.plot(sol_inhibisi.t, sol_inhibisi.y[2], color=warna['P'], lw=2.5, label='[P] Product P')
ax1.set_title('Dinamika Konsentrasi\n(dengan Allosteric Inhibition)', fontsize=11, fontweight='bold')
ax1.set_xlabel('Waktu (a.u.)')
ax1.set_ylabel('Konsentrasi (a.u.)')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# Panel 2: Perbandingan [P]
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(sol_inhibisi.t, sol_inhibisi.y[2], color=warna['P'],    lw=2.5, label='[P] dengan inhibisi')
ax2.plot(sol_baseline.t, sol_baseline.y[2], color=warna['base'], lw=2.5, ls='--', label='[P] tanpa inhibisi')
ax2.set_title('Perbandingan Akumulasi Product P\n(dengan vs tanpa Inhibisi)', fontsize=11, fontweight='bold')
ax2.set_xlabel('Waktu (a.u.)')
ax2.set_ylabel('Konsentrasi [P] (a.u.)')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# Panel 3: Laju v1
ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(sol_inhibisi.t, v1_inh,  color=warna['v1'],  lw=2.5, label='v₁ dengan inhibisi')
ax3.plot(sol_inhibisi.t, v1_base, color=warna['base'], lw=2.5, ls='--', label='v₁ tanpa inhibisi')
ax3.set_title('Laju Reaksi v₁ Sepanjang Waktu\n(Efek Feedback Inhibition)', fontsize=11, fontweight='bold')
ax3.set_xlabel('Waktu (a.u.)')
ax3.set_ylabel('Laju v₁ (a.u./waktu)')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

# Panel 4: Ringkasan semua variabel
ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(sol_inhibisi.t, sol_inhibisi.y[0], color=warna['A'],  lw=2, label='[A]')
ax4.plot(sol_inhibisi.t, sol_inhibisi.y[1], color=warna['B'],  lw=2, label='[B]')
ax4.plot(sol_inhibisi.t, sol_inhibisi.y[2], color=warna['P'],  lw=2, label='[P]')
ax4.plot(sol_inhibisi.t, v1_inh,            color=warna['v1'], lw=2, ls=':', label='v₁ (laju)')
ax4.set_title('Ringkasan Sistem:\nKonsentrasi & Laju v₁', fontsize=11, fontweight='bold')
ax4.set_xlabel('Waktu (a.u.)')
ax4.set_ylabel('Nilai (a.u.)')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)

fig.suptitle(
    'Simulasi Kinetik Jaringan Metabolik\n'
    'dengan Allosteric Non-Competitive Feedback Inhibition (P → v₁)\n'
    'BISB211605 Bioteknologi — Fakultas Biologi UGM',
    fontsize=13, fontweight='bold', y=1.01
)

plt.savefig('simulasi_metabolik.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n✅ Simulasi selesai. Grafik disimpan sebagai 'simulasi_metabolik.png'")
print(f"\nNilai steady state (approx. t=50):")
print(f"  [A] ≈ {sol_inhibisi.y[0, -1]:.4f}")
print(f"  [B] ≈ {sol_inhibisi.y[1, -1]:.4f}")
print(f"  [P] ≈ {sol_inhibisi.y[2, -1]:.4f}")
print(f"  v₁  ≈ {v1_inh[-1]:.4f}")
```

### Interpretasi Hasil Simulasi

Berdasarkan parameter yang diberikan, hasil simulasi menunjukkan tiga pola dinamika yang khas:

**Metabolit A** — Mengalami kenaikan konsentrasi yang cepat di awal karena v₁ bekerja penuh saat [P] masih mendekati nol. Setelah beberapa unit waktu, [A] mencapai **steady state** karena laju konsumsi oleh v₂ dan v₄ menyeimbangkan laju produksi dari v₁.

**Metabolit B** — Mengikuti pola yang serupa dengan A namun dengan laju yang lebih lambat karena bergantung pada konsentrasi A. B juga mencapai **steady state** setelah A mencapai keseimbangan.

**Product P** — Terus meningkat secara monoton namun dengan **laju yang semakin melambat** seiring waktu. Ini adalah bukti langsung dari mekanisme feedback inhibition: semakin tinggi [P], semakin besar faktor inhibisi (1 + [P]/Kᵢ), sehingga v₁ semakin terhambat dan seluruh flux dalam sistem secara bertahap berkurang. Fenomena ini merepresentasikan **regulasi homeostasis alami** yang dimiliki sel untuk mencegah overakumulasi produk yang bisa menjadi toksik atau membuang energi secara sia-sia.

**Laju v₁** — Tanpa inhibisi, v₁ bersifat konstan (karena [X] konstan). Dengan inhibisi, v₁ terus **menurun** seiring akumulasi P, menggambarkan self-limiting behavior yang merupakan mekanisme regulasi alami sel.

---

## 📚 Referensi

1. Palsson, B.Ø. (2015). *Systems Biology: Constraint-based Reconstruction and Analysis*. Cambridge University Press.
2. Strogatz, S.H. (2018). *Nonlinear Dynamics and Chaos*. CRC Press.
3. Fell, D.A. (1997). *Understanding the Control of Metabolism*. Portland Press.
4. GitHub Repository Matakuliah: https://github.com/lab-biotek-bio-ugm/S1_BISB211605_Biotechnology
5. Virtanen, P. et al. (2020). SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. *Nature Methods*, 17, 261–272.

---

*Dokumen ini dibuat sebagai panduan pengerjaan UAS Bioteknologi (BISB211605). Seluruh jawaban merupakan hasil pemikiran dan critical thinking mahasiswa sendiri dengan memanfaatkan sumber literasi yang tersedia sesuai ketentuan ujian.*
