# Optimizing LSPR Biosensors on Substrates

**📄 [View the poster (PDF)](poster/Optimizing_LSPR_Biosensors_on_Substrates.pdf)**

**Daniel Lee, SarahPearl Anderson, Keisuke Hasegawa**
Grinnell College — Mentored Advanced Project (MAP)

## Summary

Localized surface plasmon resonance (LSPR) occurs when the electric field of an
incident light wave excites free electrons of a noble metal nanoparticle (NP) to
collectively oscillate at a resonant wavelength. When molecules bind to an NP,
the local refractive index shifts and so does the resonant wavelength, which
makes LSPR useful for biosensors. LSPR biosensors can be used in suspension or
on a substrate; embedding NPs on substrates decreases sensitivity but is more
practical for handling, storage, and washing.

**Objective:** develop a theory to predict the sensitivity of NPs on a substrate
based on their physical dimensions.

**Methods.** All data were obtained with the MATLAB toolbox
[NANOBEM](#dependency-nanobem), an open-source Maxwell's equation solver for
nanophotonics. We tested rods (NR), triangles (NT), cylinders (NC), and disks
(ND). Discretization settings were determined by verifying data from Miller
2005, Rabbani 2024, and other works.

**Nanoparticles on substrate (one layer).** Papers such as Martinsson 2014 have
concluded that increased fractional contact area with the substrate (Xa)
linearly decreases sensitivity; our results show this is not necessarily true.
For Au NDs and NTs on glass substrate submerged in water, Xa and sensitivity
have a positive, nearly linear relationship, and for the same NPs in suspension,
aspect ratio (R) and sensitivity have a positive, hyperbolic relationship. We
concluded that higher contact area increases sensitivity for NDs and NTs because
R has a greater influence on sensitivity than Xa — indicating that flatter NP
shapes on substrate will yield more sensitive biosensors.

**Nanoparticles with an analyte layer on substrate (two layer).** We studied the
stability of NANOBEM for a silver disk with an added analyte layer on a glass
substrate, lowering the particle toward the substrate–analyte boundary to find
the minimum height at which NANOBEM gives stable sensitivities.

**Theory.** We theorize that NPs of any shape on a substrate can be modeled as a
partial sphere on the substrate, and developed an effective medium theory
relating refractive index to energy density for this geometry.

See the [poster](poster/Optimizing_LSPR_Biosensors_on_Substrates.pdf) for
figures, equations, and full conclusions.

## Code

MATLAB scripts and functions used for the NANOBEM simulations and analysis.

```
code/
├── summer-2025/
│   ├── plotPosterFigures.m          # plots the poster's S vs Xa and S vs R figures
│   ├── fminsearch-peak-finder/      # peak (resonant) wavelength search with fminsearch
│   ├── nanocylinder/  nanodisk/  nanorod/  nanotriangle/
│   │                                # per-shape simulations, spectra, sensitivity, discretization
│   ├── experimental-geometries/
│   ├── miller-2005-figures/         # reproductions used to verify discretization settings
│   ├── rabbani-2024-figures/
│   └── martinsson-2014-figures/
└── fall-2025/                       # later versions of the substrate peak-finding scripts
```

Notes:
- `summer-2025/` and `fall-2025/` contain some files with the same name
  (e.g. `subDiskLambdaMax.m`, `triangleLambdaMax.m`). Add only one of the two
  folders to the MATLAB path at a time to avoid name shadowing.
- `fall-2025/tridisklambda.m` reads a `Simulation Data.xlsx` workbook that is
  not included in this repository.

### Dependency: NANOBEM

The simulations require the **NANOBEM** toolbox by Ulrich Hohenester, which is
**not** included here. Get it from the official source and add it to your
MATLAB path:

- GitHub: <https://github.com/uhohenester/nanobem22>
- MATLAB File Exchange: <https://www.mathworks.com/matlabcentral/fileexchange/107599-nanobem22>

```matlab
addpath(genpath('/path/to/nanobem'));
addpath(genpath('code/summer-2025'));
```

Please cite NANOBEM as requested by its author (see reference 5 below).

## Acknowledgments

This work was supported by the Mentored Advanced Project (MAP) program at
Grinnell College. We thank our Professor, Keisuke Hasegawa, for his supervision
and guidance.

## References

1. Miller, Molly M., and Anne A. Lazarides. 2005. "Sensitivity of Metal Nanoparticle Surface Plasmon Resonance to the Dielectric Environment." *The Journal of Physical Chemistry B* 109 (46): 21556–65. <https://doi.org/10.1021/jp054227y>
2. Rabbani, Allison, Rachel Rudacille, and Keisuke Hasegawa. 2024. "Local Refractive Index Sensitivity of Localized Surface Plasmon Resonance Biosensors." *The Journal of Physical Chemistry C* 128 (45): 19210–21. <https://doi.org/10.1021/acs.jpcc.4c03456>
3. Martinsson, Erik, Marinus A. Otte, Mohammad Mehdi Shahjamali, Borja Sepulveda, and Daniel Aili. 2014. "Substrate Effect on the Refractive Index Sensitivity of Silver Nanoparticles." *The Journal of Physical Chemistry C* 118 (42): 24680–87. <https://doi.org/10.1021/jp5084086>
4. Mayer, Kathryn M., and Jason H. Hafner. 2011. "Localized Surface Plasmon Resonance Sensors." *Chemical Reviews* 111 (6): 3828–57. <https://doi.org/10.1021/cr100313v>
5. Hohenester, Ulrich. 2024. "Nanophotonic Resonators in Stratified Media with the Nanobem Toolbox." *Computer Physics Communications* 294 (January): 108949. <https://doi.org/10.1016/j.cpc.2023.108949>
6. Liang, Jie, Haizhou Liu, Jianyu Yu, Lin Zhou, and Jia Zhu. 2019. "Plasmon-Enhanced Solar Vapor Generation." *Nanophotonics* 8 (5): 771–86. <https://doi.org/10.1515/nanoph-2019-0039>

## License

- **Code** (`code/`): [MIT License](LICENSE)
- **Poster** (`poster/`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](poster/LICENSE)
