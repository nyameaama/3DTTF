"""
DNS-style PyFR simulation script for 3D flat-plate wing with virtual rotation.
Uses an existing wing_domain.msh mesh file (no Gmsh mesh generation).
"""

import subprocess
import math
import time
import csv
import os
import platform, shutil

# 0) Backend selection (unchanged)
backend = os.getenv('PYFR_BACKEND')
if backend is None:
    if platform.system() == 'Windows' and shutil.which('nvcc'):
        backend = 'cuda'
        precision = 'single'
    elif platform.system() == 'Darwin':
        backend = 'metal'
        precision = 'single'
    else:
        backend = 'openmp'
        precision = 'double'

print(f"Using backend: {backend}")
os.environ["OMP_NUM_THREADS"] = "10"

# 1) Mesh import (once)
mesh_file  = "wing_mesh.msh"      # ← your precomputed mesh
pyfrm_file = "wing_mesh.pyfrm"
# Convert only if .pyfrm doesn't already exist
if not os.path.exists(pyfrm_file):
    subprocess.run(
      ["pyfr", "import", mesh_file, pyfrm_file],
      check=True
    )                               # pyfr import mesh.msh mesh.pyfrm  

# 2) Simulation parameters
chord = 0.04; span = 0.40427
half_span = span / 2.0
thickness = 0.002; half_thickness = thickness/2.0

# Build case list
sweep_angles = list(range(-30,31,2))
aoa_angles   = list(range(-30,31,2))
rho_inf = 1.0; p_inf = 1.0; gamma = 1.4
Ma = 0.1; a_inf = math.sqrt(gamma*p_inf/rho_inf)
U_inf = Ma * a_inf
tot_press = p_inf * (1 + (gamma-1)/2*Ma**2)**(gamma/(gamma-1))
tot_temp  = 1 + (gamma-1)/2*Ma**2
cpTt = gamma/(gamma-1) * tot_temp

# Prepare CSV
results_csv = "aero_results.csv"
with open(results_csv, 'w', newline='') as f:
    csv.writer(f).writerow(["sweep_deg","AoA_deg","CL","CD","Cm","time_s"])

# 3) Main loop: run PyFR on the existing .pyfrm
total = len(sweep_angles)*len(aoa_angles)
count = 0
start_all = time.time()

for sweep in sweep_angles:
    for aoa in aoa_angles:
        count += 1
        # Flow-direction decomposition
        alpha = math.radians(aoa)
        beta  = math.radians(sweep)
        u_val = U_inf*math.cos(alpha)*math.cos(beta)
        v_val = U_inf*math.cos(alpha)*math.sin(beta)
        w_val = -U_inf*math.sin(alpha)

        # Write PyFR config
        cfg = f"util/case_s{int(sweep)}_a{int(aoa)}.ini"
        with open(cfg, 'w') as f:
            # ----------------------------
            # Backend & Precision
            # ----------------------------
            f.write('[backend]\n')
            f.write('precision = double\n')
            f.write('rank-allocator = linear\n\n')

            # ----------------------------
            # Physical Constants
            # ----------------------------
            f.write('[constants]\n')
            f.write('gamma = 1.4\n')
            f.write('mu = 0.001\n')
            f.write('Pr = 0.71\n\n')

            # ----------------------------
            # Solver Setup
            # ----------------------------
            f.write('[solver]\n')
            f.write('system = navier-stokes\n')
            f.write('order = 5\n\n')

            # ----------------------------
            # Time Integrator
            # ----------------------------
            f.write('[solver-time-integrator]\n')
            f.write('scheme = rk45\n')
            f.write('controller = pi\n')
            f.write('atol = 1e-4\n')
            f.write('rtol = 1e-4\n')
            f.write('tstart = 0\n')
            f.write('tend = 2\n')
            f.write('dt = 0.001\n\n')

            # ----------------------------
            # Interface & Element Settings
            # ----------------------------
            f.write('[solver-interfaces]\n')
            f.write('riemann-solver = rusanov\n')
            f.write('ldg-beta = 0.5\n')
            f.write('ldg-tau = 0.1\n\n')

            f.write('[solver-interfaces-quad]\n')
            f.write('flux-pts = gauss-legendre\n\n')

            f.write('[solver-elements-hex]\n')
            f.write('soln-pts = gauss-legendre\n\n')

            # ----------------------------
            # Triangular Interfaces
            # ----------------------------
            f.write('[solver-interfaces-tri]\n')
            f.write('flux-pts = williams-shunn\n')
            f.write('quad-deg = 1\n')
            f.write('quad-pts = williams-shunn\n\n')

            # ----------------------------
            # Tetrahedral Elements
            # ----------------------------
            f.write('[solver-elements-tet]\n')
            f.write('soln-pts = shunn-ham\n')
            f.write('quad-deg = 2\n')
            f.write('quad-pts = shunn-ham\n\n')

            # ----------------------------
            # Plugins
            # ----------------------------
            f.write('[soln-plugin-nancheck]\n')
            f.write('nsteps = 10\n\n')

            f.write('[soln-plugin-fluidforce-airfoil]\n')
            f.write('nsteps = 1\n')
            f.write('quad-deg = 4\n')
            f.write(f'file = {os.path.abspath("util/airfoil-forces.csv")}\n')
            f.write('header = true\n\n')

            # ----------------------------
            # Boundary Conditions: Inlet
            # ----------------------------
            f.write('[soln-bcs-inlet]\n')
            f.write('type = char-riem-inv\n')
            f.write(f'rho = {rho_inf}\n')
            f.write(f'u = {u_val:.12f}\n')
            f.write(f'v = {v_val:.12f}\n')
            f.write(f'w = {w_val:.12f}\n')
            f.write('p = 1\n\n')

            # ----------------------------
            # Boundary Conditions: Outlet
            # ----------------------------
            f.write('[soln-bcs-outlet]\n')
            f.write('type = sub-out-fp\n')
            f.write('p = 1\n\n')

            # ----------------------------
            # Boundary Conditions: Horizontal
            # ----------------------------
            f.write('[soln-bcs-horizontal]\n')
            f.write('type = slp-adia-wall\n\n')

            # ----------------------------
            # Boundary Conditions: Airfoil
            # ----------------------------
            f.write('[soln-bcs-airfoil]\n')
            f.write('type = no-slp-adia-wall\n\n')

            # ----------------------------
            # Initial Conditions
            # ----------------------------
            f.write('[soln-ics]\n')
            f.write('rho = 1\n')
            f.write(f'u = {u_val:.12f}\n')
            f.write(f'v = {v_val:.12f}\n')
            f.write(f'w = {w_val:.12f}\n')
            f.write('p = 1\n')

        print(f"Starting Solve:")
        # Run PyFR solve (serial or parallel)
        t0 = time.time()
        subprocess.run(
            ["pyfr", "run",
            "-b", backend,      # <-- specify the backend here
            pyfrm_file, cfg],
            check=True
        )
        dt = time.time() - t0

        # Read last line of wing_forces.csv (as before) and compute CL, CD, Cm
        with open("util/airfoil-forces.csv") as ff:
            row = list(csv.reader(ff))[-1]
            fx, fy, fz, my = map(float, (row[1], row[2], row[3], row[5]))
        # Project forces into drag/lift and compute coefficients
        V = math.sqrt(u_val**2+v_val**2+w_val**2)
        ex,ey,ez = u_val/V, v_val/V, w_val/V
        drag  = - (fx*ex + fy*ey + fz*ez)
        lift  = fz - (fx*ex + fy*ey + fz*ez)*ez
        q     = 0.5 * rho_inf * U_inf**2
        S     = chord*span
        CL    = lift / (q*S)
        CD    = drag / (q*S)
        Cm    = my   / (q*S*chord)

        # Append results
        with open(results_csv,'a',newline='') as ff:
            csv.writer(ff).writerow([sweep,aoa,CL,CD,Cm,f"{dt:.2f}"])

        # Progress
        elapsed = time.time() - start_all
        avg     = elapsed / count
        eta     = avg * (total - count)
        print(f"Case {count}/{total} (s={sweep},a={aoa}): CL={CL:.3f} CD={CD:.3f} Cm={Cm:.3f}, time={dt:.1f}s, ETA={eta:.1f}s")

print("All done.")