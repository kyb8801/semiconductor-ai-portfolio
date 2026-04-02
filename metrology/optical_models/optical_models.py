import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

lam_nm = np.linspace(200, 900, 1000)
E_eV = 1239.8 / lam_nm

def eps_to_nk(eps_complex):
    eps1 = eps_complex.real
    eps2 = eps_complex.imag
    eps_abs = np.abs(eps_complex)
    n = np.sqrt((eps_abs + eps1) / 2)
    k = np.sqrt((eps_abs - eps1) / 2)
    return n, k

def lorentz_eps(E, oscillators):
    eps = np.ones(len(E), dtype=complex)
    for (A, E0, G) in oscillators:
        eps += A / (E0**2 - E**2 - 1j*G*E)
    return eps

si_oscillators = [
    (12.0, 3.4, 0.15),
    (22.0, 4.2, 0.60),
    ( 4.0, 5.0, 1.20),
]
eps_si = lorentz_eps(E_eV, si_oscillators)
n_si, k_si = eps_to_nk(eps_si)

def tauc_lorentz_eps2(E, A, E0, Eg, G):
    eps2 = np.zeros(len(E))
    mask = E > Eg
    Em = E[mask]
    eps2[mask] = (
        A * E0 * G * (Em - Eg)**2
        / ((Em**2 - E0**2)**2 + G**2 * Em**2)
        / Em
    )
    return eps2

def tauc_lorentz_eps1_kk(E, A, E0, Eg, G):
    E_wide = np.linspace(0.01, 30.0, 8000)
    eps2_wide = tauc_lorentz_eps2(E_wide, A, E0, Eg, G)
    eps1 = np.ones(len(E))
    for i, Ei in enumerate(E):
        integrand = E_wide * eps2_wide / (E_wide**2 - Ei**2 + 1e-8)
        integrand[np.abs(E_wide - Ei) < 0.05] = 0
        eps1[i] = 1 + (2/np.pi) * np.trapz(integrand, E_wide)
    return eps1

hfo2_params = dict(A=30.0, E0=8.5, Eg=5.8, G=2.5)
eps2_hfo2 = tauc_lorentz_eps2(E_eV, **hfo2_params)
eps1_hfo2 = tauc_lorentz_eps1_kk(E_eV, **hfo2_params)
eps_hfo2  = eps1_hfo2 + 1j * eps2_hfo2
n_hfo2, k_hfo2 = eps_to_nk(eps_hfo2)

def drude_eps(E, Ep, G):
    return 1 - Ep**2 / (E**2 + 1j*G*E)

eps_al = drude_eps(E_eV, Ep=14.8, G=0.5)
n_al, k_al = eps_to_nk(eps_al)

n_sio2 = 1.458 + 0.00354 / (lam_nm/1000)**2
k_sio2 = np.zeros_like(n_sio2)
```

실행하면 `Writing optical_models.py`가 뜬다. 그 다음 cmd에서:
```
cd C:\Users\kyb88\semiconductor-ai-portfolio
copy C:\Users\kyb88\optical_models.py metrology\optical_models\
git add metrology/optical_models/optical_models.py
git commit -m "feat: add optical_models.py source code"
git push origin main
