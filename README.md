# Temperer
## Forward modeling of thermal evolution through geological time

![Build Status](https://github.com/equinor/temperer/actions/workflows/python-test.yml/badge.svg?branch=main)
![Build Status](https://github.com/equinor/temperer/actions/workflows/docs.yml/badge.svg?branch=main)

[Documentation](https://fuzzy-meme-o4w5534.pages.github.io/)

Temperer is a python package used for modeling thermal evolution based on McKenzie's type basin extension. It can be use for:

- Finding beta factor
- Calculating subsidence and thermal history
- Basement heat flow through time

## Features
- Multi-1D simulation
- Full 3D simulation with dolfinx
- Build model from either: 
    - Python objects
    - [XTGeo](https://github.com/equinor/xtgeo/) supported surface formats
- Multi-rift phase support
- Ensemble models with ERT https://github.com/equinor/ert

## Installation

Until it is available on pypi, you will need to clone the repo

```
git clone git@github.com:equinor/temperer.git
pip install .
```
For a specific release
```
git clone git@github.com:equinor/temperer.git --branch <VERSION>
pip install .
```

For full 3D simulation, dolfinx is required.

See https://docs.fenicsproject.org/dolfinx/main/python/installation.html for installation instructions.
