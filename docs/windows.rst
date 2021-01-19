Installation (Windows)
======================

Download Conda
--------------
First, a python distribution and package manager is required. For windows, we suggest |Miniconda|,
a lightweight variant of the popular Anaconda distribution. To install Miniconda, follow the |ii|.


Virtual Environment
-------------------

It is highly recommended to install LIFEsim in a *virtual environment*. This ensures installed
packages and changes made do not affect other projects on a system. The steps of creating a virtual
environment with conda are described in the following.

First, open the Conda Prompt, navigate to where you want to create the virtual environment and type

.. code-block:: console

   > conda create --name new_folder/

Some conda commands will ask for confirmation, which can be affirmed by typing.

.. code-block:: console

   Proceed ([y]/n)?
   > y

A virtual environment in the folder ``new_folder`` is created in the current directory

To activate and deactivate your new virtual environment use the following commands respectively

.. code-block:: console

   > conda activate new_folder\

.. code-block:: console

   > conda deactivate


Download from Github
--------------------

Navigate to the directory in which you desire to place the LIFEsim repository. Then clone the
repository from Github by executing

.. code-block:: console

   > git clone https://github.com/fdannert/LIFEsim.git

.. Hint:: If ``git`` is not installed on your system run

   .. code-block:: console

      > conda install git

The dependencies required by LIFEsim can be installed with

.. code-block:: console

   > conda install --file requirements_windows.txt

.. Important::
   LIFEsim need a modified version of the package |SpectRes| to run. Please install it via the
   following procedure.

   First, make sure that you are in the directory where you want to install SpectRes. Then run

   .. code-block:: console

      > git clone https://github.com/fdannert/SpectRes.git

The last step is point the Python install of your virtual environment to LIFEsim and SpectRes.
To do so, please navigate to ``site-packages`` folder of your virtual environment, most likely
located in `` C:\Users\user_name\miniconda3\envs\new_folder\Lib\site-packages``. In this directory,
create the file ``lifesim.pth`` containing the paths to LIFEsim and SpectRes separated by a new
line

.. code-block:: console

   C:\path_to_LIFEsim\LIFEsim\
   C:\path_to_SpectRes\SpectRes\


Testing the Installation
------------------------

To test the installation, open a new conda prompt and activate the virtual environment as above.
Then open Python and import LIFEsim with

.. code-block:: console

   > python

.. code-block:: python

   >>> import lifesim

If the import statement executes, the installation has been successful. As an extra test run

.. code-block:: python

   >>> lifesim.modules.constants.c
   299792000.0

This should return the speed of light in [m s
:math:`^{-1}`
].

.. Hint:: If the ``import lifesim`` command fails, the reason is likely that the the ``PYTHONPATH``
   is not set correctly. To check for this please run in Python (started with the virtual
   environment active as above)

   .. code-block:: python

      >>> import sys
      >>> sys.path

   If the path to LIFEsim ``'C:\path_to_LIFEsim\LIFEsim\'`` is not returned this is likely the
   source of the issue.

    The same test can be performed if SpectRes does not import.


.. |Miniconda| raw:: html

   <a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">Miniconda</a>

.. |ii| raw:: html

   <a href="https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/windows.html" target="_blank">installation instructions</a>

.. |SpectRes| raw:: html

   <a href="https://github.com/ACCarnall/SpectRes" target="_blank">SpectRes</a>