# glm_mda_diffusion

## or *Globule-Linker-Model, Minimum-Dissipation-Approximation diffusion coefficient calculator*

![Banner](sequence_and_shape_banner.png)

Minimum dissipation approximation is a fast algorithm for predicting the diffusive properties of intrinsically disordered proteins.

# Try with Colab

<a target="_blank" href="https://colab.research.google.com/github/RadostW/glm-mda-diffusion/blob/main/colab_notebooks/glm_mda_diffusion.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

# Installation

```bash
python3 -m pip install glm_mda_diffusion
```

# Usage as module

Basic usage:

```bash
python3 -m glm_mda_diffusion --sequence MGSS[HHHHHH]SSGLVPR
```

Sample output:
```
Computed GLM-MDA hydrodynamic radius [Ang]:
12.279165209438174
```

# Usage as package

Basic usage

```Python
import glm_mda_diffusion
glm_mda_diffusion.hydrodynamic_radius(sequence = "MGSS[HHHHHH]SSGLVPR")
```

Advanced usage (all options displayed with default values).

Options `steric_radius` and `hydrodynamic_radius` controll linker properties, while `effective_density` and `hydrdation_thickness` controll globular region properties.

```Python
import glm_mda_diffusion

glm_mda_diffusion.protein_hydrodynamic_radius(
    sequence="MGSS[HHHHHH]SSGLVPR",
    steric_radius=1.9025,  # Ang
    hydrodynamic_radius=4.2,  # Ang
    effective_density=0.52,  # Da / Ang^3
    hydration_thickness=3.0,  # Ang
    ensemble_size=30,
    bootstrap_rounds=10,
    aminoacid_masses={
        "A": 71.08,
        "C": 103.14,
        "D": 115.09,
        "E": 129.12,
        "F": 147.18,
        "G": 57.06,
        "H": 137.15,
        "I": 113.17,
        "K": 128.18,
        "L": 113.17,
        "M": 131.21,
        "N": 114.11,
        "P": 97.12,
        "Q": 128.41,
        "R": 156.2,
        "S": 87.08,
        "T": 101.11,
        "V": 99.14,
        "W": 186.21,
        "Y": 163.18,
        "Z": 0,
        "O": 0,
        "U": 0,
        "J": 0,
        "X": 0,
        "B": 0,
    },  # Da,
)
```


# License

This software is licensed under GPLv3 License

Copyright (C) Radost Waszkiewicz (2023).

# How to cite

*Hydrodynamic Radii of Intrinsically Disordered Proteins: Fast Prediction by Minimum Dissipation Approximation and Experimental Validation.*
Radost Waszkiewicz, Agnieszka Michaś, Michał K. Białobrzewski, Barbara P. Klepka, Maja K. Cieplak-Rotowska, Zuzanna Staszałek, Bogdan Cichocki, Maciej Lisicki, Piotr Szymczak, and Anna Niedźwiecka; J. Phys. Chem. Lett. (2024)

https://doi.org/10.1021/acs.jpclett.4c00312

```bibtex
@article{Waszkiewicz_2024,
	title        = {Hydrodynamic Radii of Intrinsically Disordered Proteins: Fast Prediction by Minimum Dissipation Approximation and Experimental Validation},
	author       = {Waszkiewicz, Radost and Michas, Agnieszka and Bia{\l}obrzewski, Micha{\l} K and Klepka, Barbara P and Cieplak-Rotowska, Maja K and Stasza{\l}ek, Zuzanna and Cichocki, Bogdan and Lisicki, Maciej and Szymczak, Piotr and Niedzwiecka, Anna},
	year         = 2024,
	journal      = {The Journal of Physical Chemistry Letters},
	publisher    = {ACS Publications},
	volume       = 15,
	number       = 19,
	pages        = {5024--5033}
}
```

# Bibliography

- *Diffusion coefficients of elastic macromolecules.* B. Cichocki, M. Rubin,  A. Niedzwiecka, and P. Szymczak; J. Fluid Mech. (2019)

- *GRPY: An Accurate Bead Method for Calculation of Hydrodynamic Properties of Rigid Biomacromolecules.* P. Zuk, B. Cichocki, and P. Szymczak; Biophys. J. (2018)

- *Pychastic: Precise Brownian dynamics using Taylor-Ito integrators in Python.* R. Waszkiewicz, M. Bartczak, K. Kolasa, and M. Lisicki;  SciPost Phys. Codebases (2023)
