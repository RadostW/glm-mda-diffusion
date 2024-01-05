import glm_mda_diffusion

result = glm_mda_diffusion.protein_hydrodynamic_radius(
    sequence="[GLPLPLKNENAIVDGDGTSVVTTKEDASTIFERDPNPANQVSAMVTGVILDENGD]PGESDESVENVDNDGEGGDKDDDKNGEDNDLDNKEHEEEKGDDDRGDDEEEDDAEGDNDSNDNEGDDDDDDDSGDDDDVDESGADEDDDDDSGD",
    steric_radius=1.9025,  # Ang
    hydrodynamic_radius=4.2,  # Ang
    effective_density=0.52,  # Da / Ang^3
    hydration_thickness=3.0,  # Ang
    ensemble_size=30,
    bootstrap_rounds=10,
)

print(result)
