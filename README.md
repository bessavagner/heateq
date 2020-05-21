# Heat Equation
The module heateq can be used to solve the heat equation in 1D.

## Getting Started

This module solve the partial differential (PDE) equation describing, for instance, the heat diffusion fenomena. The equation is as follows:

  <img src="https://render.githubusercontent.com/render/math?math=\LARGE \frac{\partial u}{\partial t} = D\frac{\partial^2 u}{\partial x^2}">

The boundary condition is <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\partial u}{\partial t}(-L, t) = \frac{\partial u}{\partial t}(L, t) = 0">, where <img src="https://render.githubusercontent.com/render/math?math=\large L"> is the size of the sample.

### Prerequisites

Make sure you have installed `numpy` and `scipy` in your systems or environment. It is suggested to install via pip:

```
pip install numpy
pip install scipy
```

### Usage

In example bellow, the defalt sample is initiated: 10 mm linear sample (with spacing 0.1mm - equivalent to 100 points), heated at 100°C (372 K) and smoothly cooled down until both edges at 20ºC (296 K). The choosen material is iron (D = 23 mm²/s).
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
Here's the sample image:

![Time varying map of the temperature distribution](https://github.com/bessava/heateq/blob/master/sample.png)


## Built With

* [Numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)

## Authors

* **Vagner Bessa** - *Initial work* - [bessava](https://github.com/bessava)

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/bessava/heateq/blob/master/LICENSE) file for details

## Acknowledgments

* IFCE - Instituto Federal de Ciência, Tecnologia e Educação do Ceara - Brazil 
* Inspiration: Steve Brunton [See his chanel](https://www.youtube.com/channel/UCm5mt-A4w61lknZ9lCsZtBw)
