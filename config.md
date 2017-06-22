CONFIG FILE
===========

# SOLVER SPECIFIC

## Scheme: van_leer / ldfss0 / hlle / ausm
hlle

## Higher order extension: none or ppm or muscl
ppm

## CFL
2.0

## Time-stepping method: global (g) or local (l)
g

## Higher Order Time extension: none or RK4
RK4

## Tolerance for residue norm comparisons
1e-6

## Grid file
Circle-IB-3.txt

## Immersed Boundary filename ('~' for no Immersed Boundary Method)
circle.dat

## State Load File ('~' for no load file)
###load.fvtk
~

## Max Iterations
250001

## Checkpoint iter (dump data after how many iterations)
### (Enter 0 to turn checkpointing off)
5000

## Debug level: Most detail (1) to least detail (5)
5

# FLOW SPECIFIC

## gamma (ratio of specific heats)
1.4

## R\_gas (specific gas constant)
287.

## FREE STREAM PROPERTIES

### Free Stream Density
1.225

### Free Stream X Speed
0.14

### Free Stream Y Speed
0.

### Free Stream Pressure
101325

### Viscous effects
### Give mu reference = 0 for inviscid
### Using Sutherlands law for coefficient of viscosity
### mu reference or mu0 (in kg/ms)
1.716e-5

### T reference or T0 (in K)
273.15

### Sutherland Temperature (K)
110.4

### Prandlt Number
0.7

### Post Processing

### Plot type
Pseudocolor

### Plot variable
Vorticity
