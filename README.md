
# Summary

TEMPO (Time-dependent Evolution of Multiple Pulse Operations) offers accessible and efficient simulations of pulse sequences in Python, using the suite of master equation solvers available in the Quantum Toolbox in Python (QuTiP). 
Besides straightforward definition of pulse sequence structures, including any underlying time-dependent Hamiltonians and pulse timing information, TEMPO enables faster simulations of pulse sequence dynamics (compared to naive implementations using QuTiP) while remaining compatible with the existing collection of QuTiP subpackages. Utilizing the master equation solvers that are native to QuTiP, TEMPO provides two key advantages for numerical simulations of pulse sequence dynamics.


**Ease of Use:** By incorporating both the characteristics of a Hamiltonian and its time constraints as necessary properties, pulses are constructed individually then collated to form a 'pulse sequence'. 
Time evolution is performed without the need to manually deactivate each pulse Hamiltonian outside its given time interval.
Using pulse 'recipes' in TEMPO, the creation of pulses with overlapping functional forms is streamlined, with parameters that can be tuned for individual pulses.

**Faster Executions of Time Evolution:** 
For time evolution, TEMPO organizes each pulse sequence as a series of time segments, preserving only the pulses that are active within each segment.
This avoids overheads incurred by repeated inspections of inactive pulse(s), significantly speeding up evaluation of system evolution.



## Installation


Installation can be done using package manager [pip](https://pip.pypa.io/en/stable/) and python's virtual environment module [venv](https://docs.python.org/3/library/venv.html).


1. Note the directory where **TEMPO** can be found. For this guide, it is located in ``/path/to/TEMPO``.

2. Create a virtual environment (named venv1 here), then activate the environment:
```
python -m venv venv1
source venv1/bin/activate
```
3. Install dependencies:
```
pip install -r TEMPO/requirements.txt
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

    

