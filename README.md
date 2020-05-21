# Heat Equation
The module heateq can be used to solve the heat equation in 1D.

## Getting Started

This module solve the partial differential (PDE) equation describing, for instance, the heat diffusion fenomena. The equation is as follows:

  <img src="https://render.githubusercontent.com/render/math?math=\LARGE \frac{\partial u}{\partial t} = D\frac{\partial^2 u}{\partial x^2}">

The boundary condition is <img src="https://render.githubusercontent.com/render/math?math=\large u(-L, t) = u(-L, t) = 0">, where <img src="https://render.githubusercontent.com/render/math?math=\large L"> is the size of a square sample.

### Prerequisites

Make sure you have installed `numpy` and `scipy` in your systems or environment. It is suggested to install via pip:

```
pip install numpy
pip install scipy
```

### Usage

In example bellow, the defalt sample is initiated: 10 mm square sample (with spacing 0.1mm - equivalent to 100 points), heated at 100°C (372 K) and smoothly cooled down until both edges at 20ºC (296 K). The choosen material is iron (D = 23 mm²/s).
```python
import heateq as heq
import matplotlib.pyplot as plt

D = 23  # for iron
dx = 0.1
u0, t, xscale = heq.initiate(dx=dx)
k = heq.kappa(u0, dx)  # The FFT frequencies (wavenumbers)

u = heq.solve(dudt, u0, t, D, k)

plt.figure(figsize=(8, 8))
plt.contourf(xscale, t, u, 100)
plt.colorbar()
plt.show()
```

