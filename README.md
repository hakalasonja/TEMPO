# TEMPO: A Python Package for Time Evolution of Pulse Sequences in QuTiP

<p align="center">
  <a href="https://github.com/hakalasonja/tempo/actions/workflows/test.yml">
    <img src="https://github.com/hakalasonja/tempo/actions/workflows/test.yml/badge.svg" alt="Build Status">
  </a>
</p>


# Summary

TEMPO (Time-dependent Evolution of Multiple Pulse Operations) offers accessible and efficient simulations of pulse sequences in Python, using the suite of master equation solvers available in the Quantum Toolbox in Python (QuTiP). 
Besides straightforward definition of pulse sequence structures, including any underlying time-dependent Hamiltonians and pulse timing information, TEMPO enables faster simulations of pulse sequence dynamics (compared to naive implementations using QuTiP) while remaining compatible with the existing collection of QuTiP subpackages. Utilizing the master equation solvers that are native to QuTiP, TEMPO provides two key advantages for numerical simulations of pulse sequence dynamics.


**Ease of Use:** By incorporating both the characteristics of a Hamiltonian and its time constraints as necessary properties, pulses are constructed individually then collated to form a 'pulse sequence'. 
Time evolution is performed without the need to manually deactivate each pulse Hamiltonian outside its given time interval.
Using pulse 'recipes' in TEMPO, the creation of pulses with overlapping functional forms is streamlined, with parameters that can be tuned for individual pulses.

**Faster Executions of Time Evolution:** 
For time evolution, TEMPO organizes each pulse sequence as a series of time segments, preserving only the pulses that are active within each segment.
This avoids overheads incurred by repeated inspections of inactive pulse(s), significantly speeding up evaluation of system evolution.


**Documentation** for the package can be found at [https://tempo-docs.readthedocs.io](https://tempo-docs.readthedocs.io/en/latest/).

## Installation


Installation can be done using package manager [pip](https://pip.pypa.io/en/stable/) and python's virtual environment module [venv](https://docs.python.org/3/library/venv.html).


1. Note the directory where **TEMPO** is located. For this guide, it is located in ``/path/to/TEMPO``.

2. Create a virtual environment (named venv1 here), then activate the environment:
```
python -m venv venv1
source venv1/bin/activate
```
3. Install package (pip is used here):
```
pip install /path/to/TEMPO
```
4. (Optional) If using Jupyter notebook examples found in `/path/to/TEMPO/examples`, create kernel for notebook use:
```
ipython kernel install --user --name=venv1
```


Alternatively, other environment & package managers can also be used. See below for installation steps using conda.

1. Run the following to ensure all required packages can be discovered. 
```
conda config --append channels conda-forge
```
2. Build an environment with conda, and install packages.
```
conda create --name YOUR_ENV_NAME python=3.**.**
conda activate YOUR_ENV_NAME
conda install --file /path/to/TEMPO/requirements.txt
```
Your python version should be version 3, python2 is not supported. 



## Tests

After installation, unit tests can be performed using the [pytest](https://docs.pytest.org/en/stable/) package.

1. Install pytest:
```
pip install pytest
```

2. Run tests:
```
pytest /path/to/TEMPO/tempo/tests.py

```

    

## Contributing and Support

Bug reports, questions, and feature requests are welcome through [GitHub Issues](https://github.com/hakalasonja/TEMPO/issues). For details on reporting issues, setting up a development environment, and submitting pull requests, see [CONTRIBUTING.md](CONTRIBUTING.md). For private inquiries, contact jjoon@umd.edu.



## Citing TEMPO

If you use TEMPO in your research, please cite it. Citation metadata is provided in [CITATION.cff](CITATION.cff), and GitHub's "Cite this repository" button (in the repository sidebar) generates APA and BibTeX formats from it. The BibTeX entry is:

```bibtex
@software{tempo,
  author  = {Oon, Jner Tzern and Hakala, Sonja A. and Witt, George A. and Walsworth, Ronald},
  title   = {{TEMPO}: A Python Package for Time Evolution of Pulse Sequences in {QuTiP}},
  version = {1.0.0},
  year    = {2025},
  url     = {https://github.com/hakalasonja/TEMPO}
}
```



## License

TEMPO is free software, released under the [GNU General Public License v3.0 or later](COPYING.txt) (GPL-3.0-or-later). You can redistribute and/or modify it under the terms of version 3 of the License, or (at your option) any later version published by the Free Software Foundation.