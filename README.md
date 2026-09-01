# OLB_permeability

## Purpose

`PermPorousWood3d` calculates the permeability of wood microstructures and other porous materials using the Lattice Boltzmann Method (LBM) in OpenLB 1.8.1.

The geometry is imported from a VTI (`VTK ImageData`) file generated from segmented image data. The fluid region is solved using standard BGK dynamics, while the wood material is represented as a porous medium using the Guo–Zhao porous-medium formulation with a directly specified Darcy permeability.

The code supports:

- X, Y or Z flow directions
- Parallel execution with MPI
- Direct specification of material permeability
- Automatic inlet/outlet generation
- Periodic transverse boundaries
- Velocity and pressure field export for ParaView

The implementation was developed within the EU AI-TranspWood project for permeability characterization of wood microstructures generated from image-based models.

---

# Physical Model

## Fluid Domain

Fluid voxels are simulated using:

```text
BGKdynamics
```

for materials:

```text
1 = Fluid
3 = Inlet
4 = Outlet
```

## Wood Domain

Wood voxels are simulated using:

```text
GuoZhaoBGKdynamics
```

for:

```text
5 = Porous material
```

The porous drag is controlled through a prescribed Darcy permeability:

```text
K = wallPermeability [m²]
```

rather than a porosity-based resistance model.

## Diagnostic Mode

For verification purposes, the code can optionally apply the Guo–Zhao dynamics to the whole computational domain:

```text
uniformGuoZhao = 1
```

Default:

```text
uniformGuoZhao = 0
```

which uses:

- BGK in fluid
- Guo–Zhao only inside porous material

---

# Build

## Dependencies

- OpenLB 1.8.1
- C++17 compiler
- MPI (optional)

## Example Compilation

```bash
mpicxx -std=c++17 -O3 \
    -o PermPorousWood3d \
    PermPorousWood3d.cpp \
    -lopenlb
```

Alternatively, integrate the source into the OpenLB examples build system.

---

# Command-Line Arguments

The program accepts the following arguments:

```text
PermPorousWood3d
    <vtiFile>
    <arrayName>
    <scalingFactor>
    <uPhys>
    <resolution>
    <InletPressure>
    <tau>
    <wallPermeability>
    <kinematicViscosity>
    <fluidDensity>
    <criterium>
    <flowDirection>
    [uniformGuoZhao]
```

## Parameters

### vtiFile

Input geometry file.

Example:

```text
wood_structure.vti
```

### arrayName

Name of the scalar array defining fluid voxels inside the VTI file.

Example:

```text
ImageFile
```

### scalingFactor

Physical voxel scaling factor.

Units:

```text
m
```

### uPhys

Outlet velocity magnitude.

Units:

```text
m/s
```

Must be:

```text
uPhys > 0
```

### resolution

Lattice resolution used by OpenLB.

Units:

```text
lattice cells
```

### InletPressure

Applied pressure gradient.

Units:

```text
Pa/m
```

### tau

LBM relaxation time.

Recommended:

```text
0.6 < tau < 1.2
```

### wallPermeability

Physical Darcy permeability assigned to the porous material.

Units:

```text
m²
```

Typical wood values:

```text
1e-18 ... 1e-14 m²
```

### kinematicViscosity

Fluid kinematic viscosity.

Units:

```text
m²/s
```

Example:

```text
1e-6
```

for water.

### fluidDensity

Fluid density.

Units:

```text
kg/m³
```

Example:

```text
1000
```

for water.

### criterium

Convergence tolerance used during permeability monitoring.

### flowDirection

Direction of imposed flow:

```text
0 = X
1 = Y
2 = Z
```

### uniformGuoZhao (optional)

```text
0 = Production mode
1 = Diagnostic mode
```

Default:

```text
0
```

---

# Example Run

```bash
mpirun -np 16 ./PermPorousWood3d \
    wood.vti \
    ImageFile \
    1e-5 \
    1e-3 \
    200 \
    1e5 \
    0.8 \
    1e-16 \
    1e-6 \
    1000 \
    1e-5 \
    2 \
    0
```

This example:

- loads geometry from `wood.vti`
- scales voxels by `1e-5 m`
- simulates water
- applies flow in Z direction
- uses a wood permeability of `1×10⁻¹⁶ m²`
- uses standard split dynamics (BGK + Guo–Zhao)

---

# Material Definitions

| Material | Description |
|-----------|-------------|
| 0 | Exterior |
| 1 | Fluid |
| 2 | Walls |
| 3 | Inlet |
| 4 | Outlet |
| 5 | Porous material / wood |

---

# Boundary Conditions

## Inlet

Material:

```text
3
```

Boundary type:

```text
LocalPressure
```

The inlet pressure is gradually ramped during startup to avoid numerical shocks.

## Outlet

Material:

```text
4
```

Boundary type:

```text
LocalVelocity
```

The outlet velocity magnitude is defined by:

```text
uPhys
```

## Lateral Boundaries

By default the directions orthogonal to the flow are periodic.

For example:

- flow in X → periodic Y and Z
- flow in Y → periodic X and Z
- flow in Z → periodic X and Y

---

# Geometry Processing

The VTI geometry is automatically processed as follows:

1. Read voxel data from VTI.
2. Apply physical scaling factor.
3. Add inlet/outlet padding along the flow direction.
4. Create fluid and porous-material regions.
5. Generate inlet and outlet cuboids.
6. Construct parallel domain decomposition.

---

# Simulation Time

The characteristic flow-through time is estimated as:

```text
tFlow = L / uPhys
```

where:

- `L` = domain length in flow direction

The simulation horizon is automatically set to:

```text
maxPhysT = 2 × tFlow
```

The inlet pressure ramp is applied during the initial:

```text
0.5 %
```

of the simulation time.

---

# Output

The solver generates:

- VTK/VTM visualization files
- Velocity field
- Pressure field
- Runtime statistics
- Permeability diagnostics

The generated `.vtm` master file can be loaded directly into ParaView.

---

# Visualization

Open the generated VTM file in ParaView and inspect:

- pressure field
- velocity magnitude
- velocity vectors
- flow pathways
- permeability anisotropy

Typical filters:

- Slice
- Stream Tracer
- Threshold
- Glyph

---

# Practical Recommendations

## Relaxation Time

Use:

```text
0.6 < tau < 1.2
```

for stable simulations.

## Permeability

Select permeability values representative of the porous material being studied.

Examples:

```text
Dense hardwoods     ~ 1e-18 to 1e-17 m²
Softwoods           ~ 1e-17 to 1e-15 m²
Highly permeable    > 1e-15 m²
```

## Velocity

Keep:

```text
uPhys
```

sufficiently low to remain in the Darcy regime and avoid inertial effects.

## Resolution

Increase lattice resolution when:

- velocity gradients become steep
- pore structures are small
- permeability sensitivity is high

Higher resolution generally improves permeability accuracy at increased computational cost.

---

# References

- OpenLB User Guide
- Guo, Z. and Zhao, T.S., "Lattice Boltzmann model for incompressible flows through porous media"
- AI-TranspWood project
