from setuptools import setup

with open("README.md") as f:
    long_description = f.read()

setup(
    name="glm_mda_diffusion",
    version="1.0",
    description="Minimum dissipation approximation: A fast algorithm for the prediction of diffusive properties of intrinsically disordered proteins",
    url="https://github.com/RadostW/glm-mda-diffusion/",
    author="Radost Waszkiewicz",
    author_email="radost.waszkiewicz@gmail.com",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "Documentation": "https://github.com/RadostW/glm_mda_diffusion",
        "Source": "https://github.com/RadostW/glm_mda_diffusion",
    },
    license="GNU GPLv3",
    packages=["glm_mda_diffusion"],
    install_requires=[
        "numpy",
        "scipy",
        "pygrpy",
        "sarw_spheres",
    ],
    zip_safe=False,
)
