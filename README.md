# Numerical-Ray-Tracing-of-a-Thin-Accretion-Disk-Around-a-Schwarzschild-Black-Hole

This code is released under the **MIT License**. See [`LICENSE`](LICENSE).

## Overview

This project numerically ray-traces null geodesics in the Schwarzschild spacetime to recreate the appearance of a thin, optically thick accretion disk (in the spirit of Luminet 1979). The renderer traces rays **backwards from a distant observer** and records intersections with the equatorial disk.

Key controls you’ll likely tweak first:

<img width="665" height="645" alt="Theta" src="https://github.com/user-attachments/assets/ffb6fd95-6429-4e84-a902-1732416a681b" />
- `theta_deg` — camera polar angle (tilt / elevation)
- `phi_deg` — camera azimuth angle (rotation around the black hole)
- `nx`, `ny` — output image resolution (width × height)

---

## Building

### macOS (gfortran)

1) Install a Fortran compiler  
I recommend **GNU Fortran (gfortran)**.
