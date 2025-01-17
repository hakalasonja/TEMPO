.. TEMPO documentation master file, created by
   sphinx-quickstart on Tue Dec 10 12:19:30 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TEMPO documentation
=====================================

Welcome to the TEMPO documentation page! 

Introduction
==================

TEMPO (Time-dependent Evolution of Multiple Pulse Operations) offers accessible and efficient simulations of pulse sequences in Python, using the suite of master equation solvers available in the Quantum Toolbox in Python (QuTiP). 
It enables straightforward definition of pulse sequence structures, including any underlying time-dependent Hamiltonians and pulse timing information, and faster simulations of pulse sequence dynamics (compared to naive implementations using QuTiP) while remaining compatible with the existing collection of QuTiP subpackages. Given the ubiquitous use of pulse sequences throughout quantum information/computing sciences, magnetic resonance studies, and quantum metrology, this work has immediate relevance to a wide array of research applications.


Installation
==================

Installation can be done with conda. First, ensure you have setup conda, then follow this tutorial.

1. Run the following to ensure all packages can be discovered. 

::

    conda config --append channels conda-forge


2. Build an environment with conda. Make sure you are in the root folder of this project. Your python version should be 
version 3, python 2 is not supported. 
Please ensure that the python version is compatible with the packages **numpy** and **qutip**, and that **multiprocessing** is supported. 

::

    conda create --name YOUR_ENV_NAME python=3.**.**
    conda activate YOUR_ENV_NAME
    pip install -r requirements.txt


The *requirements_docs.txt* file provides more specific packages to aid in testing and building when the library is deployed to github. 
While you are welcome to use this for installation be aware that the most streamlined installation process uses the standard *requirements.txt* file instead.

If you do not use conda, other installation approaches can work as long as you have the libraries in *requirements.txt* installed. 


Using the Docs
==================

The documentation page includes all primary modules of the library under the tab **TEMPO**. This includes the classes *Pulse*, *Pulse_recipe*, *Pulse_sequence*, *Hamiltonian*, and *Evolver*. Classes for exception handling and other miscellaneous operations are under **Utilities**. Please click on any of these modules for more details.
You may also use the available search feature in the top left of the page, or the bottom of this one. 

TEMPO Library 
==================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tempo
   misc

Indices and tables
==================

* :ref:`modindex`
* :ref:`genindex`
* :ref:`search`
