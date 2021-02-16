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

Notice the three tabs for setting the instrument options, importing and previewing the spectrum and
displaying the results.
In the *Settings* tab, the required parameters for the target star and planet can be set manually.
The parameters for the instrument can either be set manually or by pressing on a scenario button on
the left side (e.g. *Optimistic*). This then automatically sets the instruments parameters to
correspond with what the LIFE team currently views as an optimistic, baseline or pessimistic
scenario.

For importing a spectrum, navigate to the *Preview* tab.

.. image:: _static/GUI_example_spectrum.png
   :width: 70%
   :align: center

Begin by choosing a spectrum to import (in .txt format) by clicking on *Browse...*. To complete
this example, please open the file ``path_to_LIFEsim/LIFEsim/docs/_static/example_spectrum.txt`` in
the pop-up dialog. Leave the option as *absolute* to only use the imported spectrum. Setting the
option to *additive* will add the imported spectrum to the planets black body spectrum calculated
according to given parameters.

.. Hint::

   A pure black body planet can be simulated by choosing the *additive* option and leaving the
   file dialog empty.

   .. image:: _static/GUI_example_blackbody.png
      :width: 70%
      :align: center

To specify the units of the spectrum you are importing, enter them in the fields *x-axis units* and
*y-axis units*. For this example, please set *x-axis units* to ``micron`` and *y-axis units* to
``photon micron-1 s-1 m-2``.

In the *Spectrum Parameter* field the parameters used during the creation of the spectrum need to
be given. Again, for the example please set *Distance* to ``10pc``, *Planet Radius* to ``1 Earth
Radius`` and leave *Integration Time* at 0.

Pressing *Preview Spectrum* will now show the spectrum in the units specified by the user. This can
be used to check the correct import of the spectrum.

.. image:: _static/GUI_example_import.png
   :width: 70%
   :align: center

Changing the drop-down menu to *converted units* will show the spectrum in the units used in
LIFEsim. This completes setting up the simulator for a run.

Change to the *Results* tab and press *Run Simulation* on the very left. This will run the
simulation and the display the results as shown below.

.. image:: _static/GUI_example_result.png
   :width: 70%
   :align: center

Above the *Run Simulation* button you can change the *Integration Time* of the simulation and
select or deselect the inclusion of specific noise sources in the simulation.

At the bottom of the *Results* tab you can choose a location to save the results at by clicking on
*Browse...* and then save the results by clicking *Save*.

