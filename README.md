<!-- PROJECT LOGO -->
<p align="center">
  <pre>
______________________   _____ __________________   
\__    ___/\_   _____/  /     \\______   \_____  \  
  |    |    |    __)_  /  \ /  \|     ___//   |   \ 
  |    |    |        \/    Y    \    |   /    |    \
  |____|   /_______  /\____|__  /____|   \_______  /
                   \/         \/                 \/                                   
  </pre>

  <p align="center">
    Simulating overlapping pulses in a variety of shapes is a challenging task that requires repeated implementation 
    in qutip. Further, existing qutip simulations for such tasks suffer from inefficient applications of mesolve 
    with overlapping pulses. TEMPO removes these issues as a one-stop wrapper over top of qutip that streamlines and speeds up the simulation of pulse sequences. 
    <br />
    <a href="https://tempo-documentation.readthedocs.io/en/latest/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/username/tempo">View Demo</a>
    ·
    <a href="https://github.com/username/tempo/issues">Report Bug</a>
    ·
    <a href="https://github.com/username/tempo/issues">Request Feature</a>
  </p>
</p>

<!-- Badges -->
<p align="center">
  <a href="https://github.com/georgew79/tempo/actions/workflows/test.yml">
    <img src="https://github.com/georgew79/tempo/actions/workflows/test.yml/badge.svg" alt="Build Status">
  </a>
  <a href="https://github.com/username/tempo/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/username/tempo.svg" alt="License">
  </a>
  <a href="https://github.com/username/tempo/releases">
    <img src="https://img.shields.io/github/v/release/username/tempo.svg" alt="Release">
  </a>
  <a href="https://github.com/username/tempo/stargazers">
    <img src="https://img.shields.io/github/stars/username/tempo.svg" alt="Stars">
  </a>
</p>

<!-- Table of Contents -->
## Table of Contents

- [About TEMPO](#about-tempo)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## About TEMPO

**REPEATED TEXT:** Simulating overlapping pulses in a variety of shapes is a challenging task that requires repeated implementation in qutip. Further, existing qutip simulations for such tasks suffer from inefficient applications of mesolve with overlapping pulses. TEMPO removes these issues as a one-stop wrapper over top of qutip that streamlines and speeds up the simulation of pulse sequences. 

**JJ THIS SECTION IS BEST LEFT TO YOU**

## Features

- **Feature One**: Define an architecture for your pulses with 'Pulse Recipes' that allow different pulses in different scenarios to follow the same overall structure.
- **Feature Two**: Instantiate your recipes into a single sequence with our easy-to-use Hamiltonian class to define the static Hamiltonians.
- **Feature Three**: Efficiently evolve your sequence using qutip underlying code in a more efficient manner than a standard qutip simulation.

## Installation

Installation can be done quickly and easily with conda. First, ensure you have setup conda, then follow this tutorial.

1. Run the following to ensure all packages can be discovered. 

```
conda config --append channels conda-forge
```

2. Then build an environment with conda. Make sure you are in the root folder of this project. Your python version should be 
version 3, python2 is not supported. Any version should work assuming that the subsequent version is compatible with your installed
version of numpy and qutip, and multiprocessing is supported. 

```
conda create --name YOUR_ENV_NAME python=3.**.**
conda activate YOUR_ENV_NAME
pip install -r requirements.txt
```

