import numpy as np
from math import pi
from scipy.optimize import minimize
import tkinter as tk
from tkinter import simpledialog

# ==========================================
# Hilfsfunktion: komplexe Zahlen mit "i"
# ==========================================
def fmt(z, digits=6):
    re = f"{z.real:.{digits}f}"
    im = f"{abs(z.imag):.{digits}f}"
    sign = "+" if z.imag >= 0 else "-"
    return f"{re} {sign} {im} i"


# ==========================================
# Namen abfragen und Spinor φ konstruieren
# φ_k = x_k * exp(i y_k)
# x_k: Alphabetposition des k-ten Buchstabens des Vornamens
# y_k: Alphabetposition des k-ten Buchstabens des Nachnamens
# ==========================================

def letter_pos(ch: str) -> int:
    """Alphabetposition A/a=1, ..., Z/z=26. Nicht-Buchstaben -> 0."""
    ch = ch.lower()
    if 'a' <= ch <= 'z':
        return ord(ch) - ord('a') + 1
    return 0

# GUI-Root erzeugen
root = tk.Tk()
root.withdraw()  # Hauptfenster verstecken

first_name = simpledialog.askstring("Vorname", "Bitte Vornamen eingeben:")
last_name = simpledialog.askstring("Nachname", "Bitte Nachnamen eingeben:")

root.destroy()  # Tk-Fenster schließen

if not first_name or not last_name:
    raise SystemExit("Vor- und Nachname müssen eingegeben werden.")

first_name = first_name.strip()
last_name = last_name.strip()

if len(first_name) < 4 or len(last_name) < 4:
    raise SystemExit("Vor- und Nachname müssen jeweils mindestens 4 Buchstaben haben.")

# Nur die ersten 4 Buchstaben verwenden
x_vals = [letter_pos(ch) for ch in first_name[:4]]
y_vals = [letter_pos(ch) for ch in last_name[:4]]

phi = np.array([
    x_vals[k] * np.exp(1j * y_vals[k])
    for k in range(4)
], dtype=complex)

print("Erzeugter Spinor φ aus Namen:")
for i, c in enumerate(phi, 1):
    print(f"  φ[{i}] = {fmt(c)}")

phi_norm_sq = np.vdot(phi, phi).real
print("\n||φ||^2 =", phi_norm_sq)
print("x_vals (Vorname-Positionen) =", x_vals)
print("y_vals (Nachname-Positionen) =", y_vals)


# ==========================================
# Pauli-Matrizen
# ==========================================
sigma_x = np.array([[0, 1],
                    [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j],
                    [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0],
                    [0, -1]], dtype=complex)

def sigma_dot_p(p):
    px, py, pz = p
    return px * sigma_x + py * sigma_y + pz * sigma_z


# ==========================================
# u-Spinor in Standardform
#   u = sqrt(E+m) * ( χ, (σ·p/(E+m))χ )^T
# ==========================================
def u_spinor(params):
    m_raw, p_raw, th_raw, ph_raw, al_raw, be_raw = params

    m = abs(m_raw)
    p = abs(p_raw)
    theta = th_raw % pi
    phi_angle = ph_raw % (2 * pi)
    alpha = al_raw % pi
    beta = be_raw % (2 * pi)

    p_vec = np.array([
        p * np.sin(theta) * np.cos(phi_angle),
        p * np.sin(theta) * np.sin(phi_angle),
        p * np.cos(theta)
    ], dtype=float)

    E = np.sqrt(m * m + p * p)

    chi = np.array([
        np.cos(alpha / 2.0),
        np.exp(1j * beta) * np.sin(alpha / 2.0)
    ], dtype=complex)

    sqrt_E_plus_m = np.sqrt(E + m + 1e-12)

    if E + m < 1e-12:
        lower = np.zeros(2, dtype=complex)
    else:
        lower = (sigma_dot_p(p_vec) @ chi) / (E + m)

    u0 = sqrt_E_plus_m * np.concatenate([chi, lower])
    return u0, E, m, p_vec


# ==========================================
# v-Spinor in Standardform
#   v = sqrt(E+m) * ( (σ·p/(E+m))χ, χ )^T
# ==========================================
def v_spinor(params):
    m_raw, p_raw, th_raw, ph_raw, al_raw, be_raw = params

    m = abs(m_raw)
    p = abs(p_raw)
    theta = th_raw % pi
    phi_angle = ph_raw % (2 * pi)
    alpha = al_raw % pi
    beta = be_raw % (2 * pi)

    p_vec = np.array([
        p * np.sin(theta) * np.cos(phi_angle),
        p * np.sin(theta) * np.sin(phi_angle),
        p * np.cos(theta)
    ], dtype=float)

    E = np.sqrt(m * m + p * p)

    chi = np.array([
        np.cos(alpha / 2.0),
        np.exp(1j * beta) * np.sin(alpha / 2.0)
    ], dtype=complex)

    sqrt_E_plus_m = np.sqrt(E + m + 1e-12)

    if E + m < 1e-12:
        upper = np.zeros(2, dtype=complex)
    else:
        upper = (sigma_dot_p(p_vec) @ chi) / (E + m)

    v0 = sqrt_E_plus_m * np.concatenate([upper, chi])
    return v0, E, m, p_vec


# ==========================================
# Distanzfunktion mit optimaler Phase
# d(params) = min_γ 1/2 || φ - e^{iγ} ψ0(params) ||^2
# ==========================================
def distance_spinor(params, kind="u"):
    if kind == "u":
        psi0, E, m, p_vec = u_spinor(params)
    else:
        psi0, E, m, p_vec = v_spinor(params)

    norm_phi_sq = np.vdot(phi, phi).real
    norm_psi_sq = np.vdot(psi0, psi0).real

    X = np.vdot(phi, psi0)   # φ† ψ0
    abs_X = np.abs(X)

    two_d_min = norm_phi_sq + norm_psi_sq - 2 * abs_X
    return 0.5 * two_d_min


# ==========================================
# Minimierung (BFGS) mit mehreren Restarts
# ==========================================
def find_best(kind="u", restarts=20, seed=123):
    rng = np.random.default_rng(seed)
    best_res = None

    for k in range(restarts):
        guess = np.array([
            50 * rng.random(),      # m_raw
            50 * rng.random(),      # p_raw
            2 * pi * rng.random(),  # theta_raw
            2 * pi * rng.random(),  # phi_raw
            2 * pi * rng.random(),  # alpha_raw
            2 * pi * rng.random()   # beta_raw
        ])

        res = minimize(
            lambda p: distance_spinor(p, kind=kind),
            guess,
            method="BFGS",
            options={"gtol": 1e-6, "maxiter": 2000}
        )

        if best_res is None or res.fun < best_res.fun:
            best_res = res

    return best_res


print("\nStarte numerische Minimierung für u-Spinoren (BFGS) ...")
best_u = find_best(kind="u", restarts=20, seed=123)
print("u-Minimierung erfolgreich?", best_u.success)
print("Minimale Distanz d_u_min =", best_u.fun)

print("\nStarte numerische Minimierung für v-Spinoren (BFGS) ...")
best_v = find_best(kind="v", restarts=20, seed=456)
print("v-Minimierung erfolgreich?", best_v.success)
print("Minimale Distanz d_v_min =", best_v.fun)


# ==========================================
# Optimalen u- und v-Spinor + Phase berechnen
# γ_opt = -arg(φ† ψ0_opt)
# ==========================================
def build_optimal_spinor(best_res, kind="u"):
    params_opt = best_res.x
    if kind == "u":
        psi0_opt, E_opt, m_opt, p_opt = u_spinor(params_opt)
    else:
        psi0_opt, E_opt, m_opt, p_opt = v_spinor(params_opt)

    X_opt = np.vdot(phi, psi0_opt)
    gamma_opt = -np.angle(X_opt)
    psi_opt = np.exp(1j * gamma_opt) * psi0_opt

    diff = phi - psi_opt
    d_check = 0.5 * np.sum(np.abs(diff)**2)

    return psi_opt, E_opt, m_opt, p_opt, diff, d_check, params_opt


psi_u_opt, E_u, m_u, p_u, diff_u, d_u_check, params_u = build_optimal_spinor(best_u, kind="u")
psi_v_opt, E_v, m_v, p_v, diff_v, d_v_check, params_v = build_optimal_spinor(best_v, kind="v")

# Tabellen-konforme Distanzen (ohne 1/2)
d_u_table = 2 * d_u_check
d_v_table = 2 * d_v_check


# ==========================================
# Ruhe-Spin S^0 aus α,β (aus den Parametern)
# S0 = 1/2 (sinα cosβ, sinα sinβ, cosα)
# ==========================================
def rest_spin_from_params(params_opt, kind="u"):
    """
    Ruhe-Spin S^0 aus α,β:
      für u:  S^0 =  +1/2 (sinα cosβ, sinα sinβ, cosα)
      für v:  S^0 =  -1/2 (sinα cosβ, sinα sinβ, cosα)
    """
    _, _, _, _, alpha_raw, beta_raw = params_opt
    alpha = alpha_raw % np.pi
    beta  = beta_raw % (2 * np.pi)

    factor = 0.5
    if kind == "v":
        factor = -0.5

    S0x = factor * np.sin(alpha) * np.cos(beta)
    S0y = factor * np.sin(alpha) * np.sin(beta)
    S0z = factor * np.cos(alpha)
    return S0x, S0y, S0z


S0x_u, S0y_u, S0z_u = rest_spin_from_params(params_u, kind="u")
S0x_v, S0y_v, S0z_v = rest_spin_from_params(params_v, kind="v")



# ==========================================
# Ausgabe u-Spinor
# ==========================================
print("\n===== Optimaler u-Spinor =====")
for i, c in enumerate(psi_u_opt, 1):
    print(f"  u_opt[{i}] = {fmt(c)}")

print("\nPhysikalische Größen (u):")
print("  E_u =", E_u)
print("  m_u =", m_u)
print("  p_u =", p_u)
print("  ||u_opt||^2 =", np.vdot(psi_u_opt, psi_u_opt).real)
print("  Distanz d_u (Definition)       =", d_u_check)
print("  Distanz d_u (Tabellenkonvention) =", d_u_table)

print("\nRuhe-Spin S^0 (u):")
print("  Sx^0 =", S0x_u)
print("  Sy^0 =", S0y_u)
print("  Sz^0 =", S0z_u)
print("  |S^0| =", 0.5)

# ==========================================
# Ausgabe v-Spinor
# ==========================================
print("\n===== Optimaler v-Spinor =====")
for i, c in enumerate(psi_v_opt, 1):
    print(f"  v_opt[{i}] = {fmt(c)}")

print("\nPhysikalische Größen (v):")
print("  E_v =", E_v)
print("  m_v =", m_v)
print("  p_v =", p_v)
print("  ||v_opt||^2 =", np.vdot(psi_v_opt, psi_v_opt).real)
print("  Distanz d_v (Definition)       =", d_v_check)
print("  Distanz d_v (Tabellenkonvention) =", d_v_table)

print("\nRuhe-Spin S^0 (v):")
print("  Sx^0 =", S0x_v)
print("  Sy^0 =", S0y_v)
print("  Sz^0 =", S0z_v)
print("  |S^0| =", 0.5)

# ==========================================
# Entscheidung: welcher ist näher?
# ==========================================
print("\n===== Vergleich =====")
print("  d_u (Def.) =", d_u_check)
print("  d_v (Def.) =", d_v_check)

if d_u_check < d_v_check:
    print("\nErgebnis: φ ist näher an einem u-Spinor als an einem v-Spinor.")
else:
    print("\nErgebnis: φ ist näher an einem v-Spinor als an einem u-Spinor.")
