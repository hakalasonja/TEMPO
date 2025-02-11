.. TEMPO documentation master file, created by
   sphinx-quickstart on Tue Dec 10 12:19:30 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TEMPO documentation
=====================================

Welcome to the TEMPO documentation page! 


TEMPO (Time-dependent Evolution of Multiple Pulse Operations) offers accessible and efficient simulations of pulse sequences in Python, using the suite of master equation solvers available in the Quantum Toolbox in Python (QuTiP). 
It enables straightforward definition of pulse sequence structures, including any underlying time-dependent Hamiltonians and pulse timing information, and faster simulations of pulse sequence dynamics (compared to naive implementations using QuTiP) while remaining compatible with the existing collection of QuTiP subpackages. Given the ubiquitous use of pulse sequences throughout quantum information/computing sciences, magnetic resonance studies, and quantum metrology, this work has immediate relevance to a wide array of research applications.


Using the Docs
==================

The documentation page includes all primary modules of the library under the tab **TEMPO**. This includes the classes *Pulse*, *Pulse_recipe*, *Pulse_sequence*, *Hamiltonian*, and *Evolver*. Classes for exception handling and other miscellaneous operations are under **Utilities**. Please click on any of these modules for more details.
You may also use the available search feature in the top left of the page, or the bottom of this one. 

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   tutorials
   tempo
   misc

Index
==================

* :ref:`modindex`
* :ref:`genindex`
* :ref:`search`
