# ITACA Detector Design Parameters

## 1. Fiducial Volume

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Fiducial diameter | D_fid | 3.20 | m |
| Fiducial radius | R_fid | 160.0 | cm |
| Drift length (TPC) | L_TPC | 1.50 | m |
| Fiducial volume | V_fid | 12.1 | m³ |
| Fiducial mass | m_fid | 1050 | kg (1.05 ton) |

## 2. Operating Conditions

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Operating pressure | P | 15 | bar |
| Xenon density | ρ_Xe | 87.0 | kg/m³ |
| Temperature | T | 300 | K |
| Ion species | — | Xe₂⁺ | — |

## 3. Electric Field & Drift

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Drift field | E | 300 | V/cm |
| Ion drift velocity | v_ion | 15 | cm/s |
| Electron drift velocity | v_e | 1.0 | mm/μs |
| Ion reduced mobility (STP) | μ_STP | 0.74 | cm²/V/s |
| Ion drift time (full) | t_ion | 10.0 | s |
| Electron drift time (full) | t_e | 1.5 | ms |

## 4. Diffusion

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Ion transverse diffusion (L=150 cm) | σ_T,ion | ~1.6 | mm |
| Electron diffusion coefficient | σ_coeff | 3.5 | mm×√(bar/cm) |
| Electron transverse diffusion (L=150 cm, P=15 bar) | σ_T,e | ~11.1 | mm |

## 5. Voltages

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Cathode voltage | V_c | 0.5 | kV |
| FAT-GEM voltage | V_FATGEM | 15.0 | kV |
| Electrode 1 (at L=1.5m, E=300 V/cm) | V_1 | 45.5 | kV |
| Electrode 2 | V_2 | 60.5 | kV |

## 6. Cathode Mesh

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Wire diameter | w | 200 | μm |
| Mesh pitch | M | 5 | mm |
| Transparency | f | 92.2 | % |

## 7. Radial Structure (center outward)

| Component | Thickness | Cumulative R | Unit |
|-----------|-----------|--------------|------|
| TPC (fiducial) | — | 160.0 | cm |
| BFD (Buffer Field) | 1.0 | 161.0 | cm |
| Clearances | 1.0 | 162.0 | cm |
| ICS (Inner Copper Shield) | 15.0 | 177.0 | cm |
| PV wall | 1.11 | 178.1 | cm |

## 8. Axial Structure (top to bottom)

| Component | Thickness | Unit |
|-----------|-----------|------|
| ICS (top) | 15 | cm |
| DSP + Light Guide + FAT-GEM | 10 | cm |
| TPC (drift region) | 150 | cm |
| Cathode | ~0 | cm |
| MRS area | 10 | cm |
| ICS (bottom) | 15 | cm |
| **Total internal** | **200** | **cm** |

## 9. Pressure Vessel (Ti Grade 5)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Material | — | Ti-6Al-4V (Grade 5) | — |
| Allowable stress | S | 240 | MPa |
| Joint efficiency | E_w | 1.0 | — |
| Inner diameter | D_i | 3.540 | m |
| Outer diameter | D_o | 3.562 | m |
| Cylinder wall thickness | t_cyl | 11.10 | mm |
| Hemisphere wall thickness | t_hem | 5.53 | mm |
| Internal height | L_i | 2.00 | m |
| External height (with caps) | L_o | 2.01 | m |

## 10. MARS Kinematics

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Arm radius | R | 1.6 | m |
| Arm linear mass density | μ | 0.24 | kg/m |
| Ion plate mass | m_plate | 0.6 | kg |
| Motor torque | τ_motor | 140 | N·m |
| Rotation angle | Δθ | π/2 | rad |
| Gas density | ρ | 87 | kg/m³ |
| Drag coefficient | C_d | 0.1 | — |
| Wake diameter | d_wake | 9.6 | mm |
| Settling exponent | n | 1.2 | — |
| Rotation time | t_rot | ~0.42 | s |
| Settling time (v_d=0.15 m/s) | t_settle | ~0.17 | s |
| Cycle time | t_cycle | ~0.59 | s |
| Dead zone | Z_dead | ~8.9 | cm |
| Fiducial efficiency | ε_geo | ~94 | % |

## 11. MARS Arm Structural (Ti Grade 5)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Number of longerons | N | 3 | — |
| Tube outer diameter | d_o | 8 | mm |
| Tube wall thickness | t_w | 0.8 | mm |
| Longeron positions | y | [-30, 0, +30] | mm |
| Ti density | ρ_Ti | 4430 | kg/m³ |
| Young's modulus | E_Ti | 114 | GPa |
| Gravity | g | 9.81 | m/s² |
| Ti Gr.5 yield strength | σ_yield | 830 | MPa |
| Allowable stress (SF=2) | σ_allow | 400 | MPa |

## 12. Physical Constants

| Constant | Symbol | Value | Unit |
|----------|--------|-------|------|
| Boltzmann constant | k_B | 1.380649×10⁻²³ | J/K |
| Elementary charge | q_e | 1.602176634×10⁻¹⁹ | C |
