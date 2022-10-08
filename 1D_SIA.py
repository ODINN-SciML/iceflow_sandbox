""" Simple glacier model.

License: BSD3

Copright: OGGM developpers, http://oggm.org
"""

import numpy as np
import matplotlib.pyplot as plt

sec_in_day = 24 * 60 * 60
sec_in_year = sec_in_day * 365

glen_n = 3
rho = 900
g = 9.81


def glacier_evolution(
    dx=100,  # grid resolution in m
    nx=200,  # grid size
    width=600,  # glacier width in m
    top_h=3000,  # bed top altitude
    bottom_h=1200,  # bed bottom altitude
    glen_a=2.4e-24,  # ice stiffness
    ela_h=2600,  # mass balance model Equilibrium Line Altitude
    mb_grad=3,  # linear mass balance gradient (unit: [mm w.e. yr-1 m-1])
    n_years=700,  # simulation time in years
):

    bed_h = np.linspace(top_h, bottom_h, nx)
    surface_h = bed_h.copy()
    thick = bed_h * 0

    def get_mb(heights):
        mb = (heights - ela_h) * mb_grad
        return mb / sec_in_year / rho

    t = 0
    dt = sec_in_day * 10

    years = np.arange(0, n_years+1, dtype=np.int64)
    volume = np.empty(n_years + 1, dtype=np.float64)
    length = np.empty(n_years + 1, dtype=np.float64)

    for i, y in enumerate(years):
        end_t = y * sec_in_year
        # Time integration
        while t < end_t:

            # This is to guarantee a precise arrival on a specific date if asked
            remaining_t = end_t - t
            if remaining_t < dt:
                dt = remaining_t

            # Surface gradient
            surface_gradient = np.zeros(nx)
            surface_gradient[1:nx-1] = (surface_h[2:] - surface_h[:nx-2])/(2*dx)
            surface_gradient[[0, -1]] = 0  # no flux boundary conditions
            
            # Diffusivity
            diffusivity = width * (rho*g)**3 * thick**3 * surface_gradient**2
            diffusivity *= 2/(glen_n+2) * glen_a * thick**2
            #print(np.mean(diffusivity))
            # Ice flux is computed on staggered gridd
            diffusivity_s = np.zeros(nx)
            diffusivity_s[1:] = (diffusivity[:nx-1] + diffusivity[1:])/2.
            diffusivity_s[0] = diffusivity[0]
            
            surface_gradient_s = np.zeros(nx)
            surface_gradient_s[1:] = (surface_h[1:] - surface_h[:nx-1])/dx
            surface_gradient_s[0] = 0

            grad_x_diff = surface_gradient_s * diffusivity_s
            flux_div = (grad_x_diff[1:] - grad_x_diff[:nx-1]) / dx

            # Mass balance
            mb = get_mb(surface_h[:nx-1])

            # Ice thickness update: old + flux div + mb
            new_thick = np.zeros(nx)
            new_thick[:nx-1] = thick[:nx-1] + (dt / width) * flux_div + dt * mb
            dH = (1.0 / width) * flux_div + 1.0 * mb
            #print(flux_div[0])
            
            if thick[-2] > 0:
                raise RuntimeError('Glacier exceeding boundaries')

            # We can have negative thickness because of MB - correct here
            thick = np.clip(new_thick, 0, None)

            # Prepare for next step
            surface_h = bed_h + thick
            t += dt

        volume[i] = (thick * width * dx).sum()
        length[i] = (thick > 0).sum() * dx

    # xcoordinates
    xc = np.arange(nx) * dx

    return xc, bed_h, surface_h, years, volume, length

xc, bed_h, surface_h, years, volume, length = glacier_evolution(n_years=700)


import time

start = time.time()
xc, bed_h, surface_h, years, volume, length = glacier_evolution(n_years=700)

print(f'Time: {time.time() - start}')
