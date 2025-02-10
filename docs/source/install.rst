============
Installation
============

Installation can be done using package manager **pip** and python's virtual environment module **venv**.

1. Navigate to parent directory, if TEMPO is located in `/path/to/TEMPO`
.. code-block:: bash

    cd /path/to/

2. Create virtual environment named venv1 in parent directory, then activate the environment::

    python -m venv venv1
    source venv1/bin/activate

3. Install dependencies::
    
    pip install -r TEMPO/requirements.txt
    
4. (Optional) If using Jupyter notebook examples found in `/path/to/TEMPO/examples`, create kernel for notebook use: ::

    ipython kernel install --user --name=venv1
