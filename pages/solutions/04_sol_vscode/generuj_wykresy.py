from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- KONFIGURACJA DANYCH Z INDEKSU ---
INDEX_FIRST_TWO = 33  # Podmień na 2 pierwsze cyfry swojego indeksu
INDEX_LAST_TWO = 24   # Podmień na 2 ostatnie cyfry swojego indeksu

v0 = float(INDEX_FIRST_TWO)
alpha_deg = float(INDEX_LAST_TWO)
alpha_rad = np.radians(alpha_deg)
g = 9.81
m = 1.0       # masa [kg]
b = 0.15      # współczynnik oporu powietrza [kg/s]

OUT = Path("files/rzut_ukosny")
OUT.mkdir(parents=True, exist_ok=True)

# Styl wykresów spójny z przykładowym raportem
INK = "#233044"
BLUE, ORANGE, GREEN, PURPLE = "#2463A6", "#D47722", "#168578", "#8455A4"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.titlesize": 13, "axes.titleweight": "bold",
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.18,
    "lines.linewidth": 2.3, "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

def zapisz(fig, nazwa):
    fig.savefig(OUT / f"{nazwa}.png", dpi=180)
    plt.close(fig)

# --- 1. SCHEMAT FIZYCZNY ---
fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
ax.set(xlim=(-0.5, 4), ylim=(-0.5, 3))
ax.axis("off")
ax.annotate("", xy=(3.5, 0), xytext=(-0.2, 0), arrowprops=dict(arrowstyle="->", lw=1.5, color=INK))
ax.annotate("", xy=(0, 2.7), xytext=(0, -0.2), arrowprops=dict(arrowstyle="->", lw=1.5, color=INK))
ax.text(3.6, -0.1, "x (poziom)", fontsize=11)
ax.text(-0.3, 2.7, "y (pion)", fontsize=11)
ax.plot(0, 0, 'o', color=BLUE, ms=8)
ax.annotate("", xy=(1.8*np.cos(alpha_rad), 1.8*np.sin(alpha_rad)), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", lw=2.5, color=ORANGE))
ax.text(1.9*np.cos(alpha_rad), 1.9*np.sin(alpha_rad), r"$v_0$", color=ORANGE, fontsize=14, weight="bold")
arc = np.linspace(0, alpha_rad, 50)
ax.plot(0.6*np.cos(arc), 0.6*np.sin(arc), color=GREEN, lw=1.5)
ax.text(0.7*np.cos(alpha_rad/2), 0.7*np.sin(alpha_rad/2), r"$\alpha$", color=GREEN, fontsize=13)
ax.set_title("Schemat rzutu ukośnego", loc="left")
zapisz(fig, "schemat")

# --- 2. SOLVER NUMERYCZNY DLA OBU MODELI ---
# Bez oporu (analitycznie)
t_fly_bez = 2 * v0 * np.sin(alpha_rad) / g
t_bez = np.linspace(0, t_fly_bez, 500)
x_bez = v0 * np.cos(alpha_rad) * t_bez
y_bez = v0 * np.sin(alpha_rad) * t_bez - 0.5 * g * t_bez**2
vx_bez = np.full_like(t_bez, v0 * np.cos(alpha_rad))
vy_bez = v0 * np.sin(alpha_rad) - g * t_bez

# Z oporem (numerycznie)
def rhs(t, s):
    x, y, vx, vy = s
    v = np.hypot(vx, vy)
    return [vx, vy, -(b/m)*vx, -g - (b/m)*vy]

def hit_ground(t, s):
    return s[1]
hit_ground.terminal = True
hit_ground.direction = -1

sol = solve_ivp(rhs, (0, 100), [0, 0, v0*np.cos(alpha_rad), v0*np.sin(alpha_rad)],
                events=hit_ground, rtol=1e-9, atol=1e-9, max_step=0.01)

t_z = sol.t
x_z, y_z, vx_z, vy_z = sol.y

# --- 3. TRAJEKTORIA PORÓWNANIE ---
fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
ax.plot(x_bez, y_bez, color=BLUE, label="Bez oporu (idealny)")
ax.plot(x_z, y_z, color=ORANGE, ls="--", label=f"Z oporem powietrza (b = {b} kg/s)")
ax.set(xlabel="Odległość x [m]", ylabel="Wysokość y [m]", title="Trajektoria ruchu w rzucie ukośnym")
ax.legend()
zapisz(fig, "trajektoria_porownanie")

# --- 4. PRĘDKOŚĆ W FUNKCJI CZASU ---
fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
axes[0].plot(t_bez, vx_bez, color=BLUE, label="$v_x$ bez oporu")
axes[0].plot(t_z, vx_z, color=ORANGE, ls="--", label="$v_x$ z oporem")
axes[0].set_ylabel("Prędkość $v_x$ [m/s]")
axes[0].legend()

axes[1].plot(t_bez, vy_bez, color=BLUE, label="$v_y$ bez oporu")
axes[1].plot(t_z, vy_z, color=ORANGE, ls="--", label="$v_y$ z oporem")
axes[1].set(xlabel="Czas t [s]", ylabel="Prędkość $v_y$ [m/s]")
axes[1].legend()
axes[0].set_title("Składowe prędkości w czasie", loc="left")
zapisz(fig, "predkosc_czas")

# --- 5. ENERGIA W FUNKCJI CZASU ---
Ek_bez = 0.5 * m * (vx_bez**2 + vy_bez**2)
Ep_bez = m * g * y_bez
E_bez = Ek_bez + Ep_bez

Ek_z = 0.5 * m * (vx_z**2 + vy_z**2)
Ep_z = m * g * y_z
E_z = Ek_z + Ep_z

fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
ax.plot(t_bez, E_bez, color=BLUE, label="Całkowita (bez oporu)")
ax.plot(t_z, E_z, color=ORANGE, ls="--", label="Całkowita (z oporem)")
ax.plot(t_z, Ek_z, color=GREEN, ls=":", label="Kinetyczna (z oporem)")
ax.plot(t_z, Ep_z, color=PURPLE, ls=":", label="Potencjalna (z oporem)")
ax.set(xlabel="Czas t [s]", ylabel="Energia [J]", title="Bilans energii mechanicznej")
ax.legend()
zapisz(fig, "energia_czas")

# --- 6. PRZESTRZEŃ FAZOWA (y vs vy) ---
fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
ax.plot(y_bez, vy_bez, color=BLUE, label="Bez oporu")
ax.plot(y_z, vy_z, color=ORANGE, ls="--", label="Z oporem")
ax.set(xlabel="Wysokość y [m]", ylabel="Prędkość pionowa $v_y$ [m/s]", title="Portret fazowy w osi pionowej (y, v_y)")
ax.legend()
zapisz(fig, "przestrzen_fazowa")

# --- 7. WPŁYW WSPÓŁCZYNNIKA OPORU b ---
fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
ax.plot(x_bez, y_bez, color=BLUE, label="b = 0 (bez oporu)")
for b_val, color in zip([0.05, 0.15, 0.30], [GREEN, ORANGE, PURPLE]):
    sol_b = solve_ivp(lambda t, s: [s[2], s[3], -(b_val/m)*s[2], -g - (b_val/m)*s[3]],
                      (0, 100), [0, 0, v0*np.cos(alpha_rad), v0*np.sin(alpha_rad)],
                      events=hit_ground, max_step=0.01)
    ax.plot(sol_b.y[0], sol_b.y[1], color=color, ls="--", label=f"b = {b_val} kg/s")
ax.set(xlabel="Odległość x [m]", ylabel="Wysokość y [m]", title="Wpływ współczynnika oporu powietrza b na zasięg")
ax.legend()
zapisz(fig, "wplyw_b")

print("Wykresy zostały pomyślnie wygenerowane w katalogu files/rzut_ukosny/")