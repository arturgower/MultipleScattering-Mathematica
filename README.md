
# MultipleScattering2D.wl - exact multiple scattering in 2D for Mathematica

A package to calculate multiple scattering according to the 2D wave equation. The best place to get started is the notebook [TwoBodyScattering.nb](examples/TwoBodyScattering.nb), though there is some 
[documention](Readme.pdf).


At present all scatterers need to be cylinders of the same size, but can be placed anywhere. The incident wave can be anything, but the default for the package is the 2D green's function. There are functions to calculate the total wave in frequency and time, together with examples on how to plot them.

## Examples 

By lining up cylinders above and below the source, we effectively create two walls:

![TwoWalls](media/TwoWallsBodyScattering.gif)


### The package's main focus is to study scattering from randomly placed cyclinders:
#### Dirichlet boundary conditions

![dirichlet](media/45-Wave_10-Scatterers.GIF)

#### Neumann boundary conditions

![neumann](media/45-Wave_12-Scatterers-Neuman.GIF)
