First Example
=============

The following is a minimal working example of simulating an observation of an exoplanet with LIFE.

As a first step, activate the virtual environment containing the LIFEsim install with

.. code-block:: console

   $ source path_to_new_folder/new_folder/bin/activate

Then, open the LIFEsim spectrum simulator GUI by running

.. code-block:: console

   $ python path_to_LIFEsim/LIFEsim/lifesim/gui/spectrum_gui.py

With executing this command, the following GUI will open

.. image:: _static/GUI_example.png
   :width: 70%
   :align: center

Now, open the file ``path_to_LIFEsim/LIFEsim/docs/_static/example_spectrum.txt`` by clicking on
*Browse...* under *Settings -> Planet -> Spectrum*.

The imported spectrum can be previewed by pressing the *Preview Spectrum* button in the *Preview*
tab. The following will be displayed

.. image:: _static/GUI_example_spectrum.png
   :width: 70%
   :align: center

The simulation can be run by pressing the *Run Simulation* button on the very left. This will
display the results in the *Results* tab.

