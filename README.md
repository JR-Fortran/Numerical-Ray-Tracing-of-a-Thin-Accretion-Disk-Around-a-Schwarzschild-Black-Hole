# Numerical-Ray-Tracing-of-a-Thin-Accretion-Disk-Around-a-Schwarzschild-Black-Hole

This code is released under the **MIT License**. See [`LICENSE`](LICENSE).

## Overview

This project numerically ray-traces null geodesics in the Schwarzschild spacetime to recreate the appearance of a thin, optically thick accretion disk (in the spirit of Luminet 1979). The renderer traces rays **backwards from a distant observer** and records intersections with the equatorial disk.

Key controls you’ll likely tweak first:

![Changing theta_deg tilts the camera](Images/Theta.png)
- `theta_deg` — camera polar angle (tilt / elevation)
![Changing phi_deg rotates the view](Images/Phi.png)
- `phi_deg` — camera azimuth angle (rotation around the black hole)
- `nx`, `ny` — output image resolution (width × height)

---

## Building

### macOS (gfortran)

1) Install a Fortran compiler  
I recommend **GNU Fortran (gfortran)**.
