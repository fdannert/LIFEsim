Installation
============

LIFEsim is available for installation from |github|. It is compatible with Python 3.7.

Virtual Environment
-------------------
.. TODO: link pip, add installation for macos and windows

It is highly recommended to install LIFEsim in a *virtual environment*. This ensures installed
packages and changes made do not affect other projects on a system. The package dependencies of
LIFEsim are best managed with a package manager like pip.

As a first step, *virtualenv* is installed with pip

.. code-block:: console

   $ pip install virtualenv

Then a virtual environment in the folder ``new_folder`` is created in the current directory

.. code-block:: console

   $ virtualenv -p python3 new_folder

To activate and deactivate your new virtual environment use the following statement respectively

.. code-block:: console

   $ source path_to_new_folder/new_folder/bin/activate

.. code-block:: console

   $ deactivate


Download from Github
--------------------

Navigate to the directory in which you desire to place the LIFEsim repository. Then clone the
repository from Github by executing

.. code-block:: console

   $ git clone https://github.com/fdannert/LIFEsim.git

.. Hint:: If ``git`` is not installed on your system run

   .. code-block:: console

      $ sudo apt install git

The dependencies required by LIFEsim can be installed with

.. code-block:: console

   $ pip install -r LIFEsim/requirements.txt

To upgrade already installed dependencies to LIFEsim requirements run

.. code-block:: console

   $ pip install --upgrade -r LIFEsim/requirements.txt

The last step is point the Python install of your virtual environment to LIFEsim. Please do so by
running

.. code-block:: console

   $ echo "export PYTHONPATH='$PYTHONPATH:/path_to_LIFEsim/LIFEsim/'" >> path_to_new_folder/new_folder/bin/activate

Testing the Installation
------------------------

To test the installation, open a new console and activate the virtual environment as above. Then
open Python and import LIFEsim with

.. code-block:: console

   $ python

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

   If the path to LIFEsim ``'/path_to_LIFEsim/LIFEsim/'`` is not returned in the results, please
   open the file ``path_to_new_folder/new_folder/bin/activate`` with a text editor of your choice.
   Then make sure that the last line of the file reads

   .. code-block:: none

      export PYTHONPATH=':/path_to_LIFEsim/LIFEsim/'


.. |github| raw:: html

   <a href="https://github.com/fdannert/LIFEsim" target="_blank">Github</a>