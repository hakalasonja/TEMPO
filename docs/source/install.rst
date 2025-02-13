Installation
============

Installation can be done using package manager `pip <https://pip.pypa.io/en/stable/>`_ and python's virtual environment module `venv <https://docs.python.org/3/library/venv.html>`_.


1. Note the directory where **TEMPO** can be found. For this guide, it is located in ``/path/to/TEMPO``.

2. Create a virtual environment (named venv1 here), then activate the environment::

    python -m venv venv1
    source venv1/bin/activate

3. Install dependencies::
    
    pip install -r /path/to/TEMPO/requirements.txt
    
4. (Optional) If using Jupyter notebook examples found in ``/path/to/TEMPO/examples``, create kernel for notebook use::

    ipython kernel install --user --name=venv1


Tests
============

After installation, unit tests can be performed using the `pytest <https://docs.pytest.org/en/stable/>`_ package.

1. Install pytest::

    pip install pytest

2. Run tests::

    pytest /path/to/TEMPO/tempo/tests.py


