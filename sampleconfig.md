CONFIG FILE
===========

# SOLVER SPECIFIC

## Scheme: van_leer / ldfss0 / hlle / ausm
hlle

## Higher order extension: none or ppm or muscl
ppm

## CFL
0.8

## Time-stepping method: global (g) or local (l)
l

## Higher Order Time extension: none or RK4
RK4

## Tolerance for residue norm comparisons
1e-6

## Grid file
expansion_fan.txt

## Immersed Boundary filename ('~' for no Immersed Boundary Method)
~

## State Load File ('~' for no load file)
###load.fvtk
~

## Max Iterations
30000

## Checkpoint iter (dump data after how many iterations)
### (Enter 0 to turn checkpointing off)
500

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
680.588

### Free Stream Y Speed
0.

### Free Stream Pressure
101325

### Viscous effects
### Give mu reference = 0 for inviscid
### Using Sutherlands law for coefficient of viscosity
### mu reference or mu0 (in kg/ms)
0.0

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
Mach_No
