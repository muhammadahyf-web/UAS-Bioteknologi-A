import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parameter dari tabel soal
V1max = 5.0   # Vmax reaksi v1
Km1   = 2.0   # Michaelis constant v1
Ki    = 3.0   # Inhibition constant
X     = 10.0  # Konsentrasi substrat eksternal (konstan)
k2    = 1.0   # Rate constant v2: A → B
k3    = 0.8   # Rate constant v3: B → P
k4    = 0.3   # Rate constant v4: A → Byproduct

def sistem_ode(t, y):
    A, B, P = y

    # v1 dengan non-competitive inhibition oleh P
    v1 = (V1max * X) / ((Km1 + X) * (1 + P / Ki))
    v2 = k2 * A
    v3 = k3 * B
    v4 = k4 * A

    dA_dt = v1 - v2 - v4
    dB_dt = v2 - v3
    dP_dt = v3

    return [dA_dt, dB_dt, dP_dt]

# Kondisi awal: semua metabolit = 0 pada t=0
y0     = [0.0, 0.0, 0.0]
t_span = (0, 50)
t_eval = np.linspace(0, 50, 1000)

# Selesaikan sistem ODE secara numerik
sol = solve_ivp(sistem_ode, t_span, y0, t_eval=t_eval, method='RK45')

# Visualisasi hasil simulasi
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='[A] Intermediate A', color='royalblue',  linewidth=2)
plt.plot(sol.t, sol.y[1], label='[B] Intermediate B', color='darkorange', linewidth=2)
plt.plot(sol.t, sol.y[2], label='[P] Product P',      color='seagreen',   linewidth=2)

plt.xlabel('Waktu (a.u.)', fontsize=12)
plt.ylabel('Konsentrasi (a.u.)', fontsize=12)
plt.title('Simulasi Kinetik Jaringan Metabolik\ndengan Allosteric Feedback Inhibition (P → v1)', fontsize=13)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('simulasi_metabolik.png', dpi=150)
plt.show()
