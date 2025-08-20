"""
TEMPO: A Python Package for Time Evolution of Pulse Sequences in QuTiP

TEMPO (Time-dependent Evolution of Multiple Pulse Operations) offers accessible and 
efficient simulations of pulse sequences in Python, using the suite of master equation 
solvers available in the Quantum Toolbox in Python (QuTiP).
"""

__version__ = "1.0.0"
__author__ = "Sonja Hakala, Jner Tzern Oon, George A. Witt, Ronald Walsworth"
__email__ = "hakalas@terpmail.umd.edu"

# Import main classes and functions that users need
from .hamiltonian import Hamiltonian
from .pulse_recipe import Pulse_recipe
from .pulse import Pulse
from .pulse_sequence import Pulse_sequence
from .evolver import Evolver


# Import any exceptions
from .exceptions import *

# Define what gets imported with "from tempo import *"
__all__ = [
    'Hamiltonian',
    'Pulse_recipe', 
    'Pulse',
    'Pulse_sequence',
    'Evolver',
    '__version__'
    ]

