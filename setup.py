from setuptools import setup

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="astro-tigress",
    version="1.0",
    description="",
    author="",
    author_email="",
    packages=["astro_tigress"],
    install_requires=install_requires,
    package_data={"astro_tigress": ["coolftn.txt"]},
)
