from setuptools import setup

setup(
    name="astro-tigress",
    version="1.0",
    description="",
    author="",
    author_email="",
    packages=["astro_tigress"],
    install_requires=["jupyter-book", "matplotlib", "numpy", "scipy", "astropy",
                      "pandas", "docutils==0.17.1", "yt", "xarray", "h5py", "urllib3",
                      "tqdm"],
    package_data={"astro_tigress": ["coolftn.txt"]},
)
