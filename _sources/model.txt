Model Processing
=================

.. autoclass:: pyfuse.pyFUSE
   :members: get_all_runids, clean_outputs, load_new_pars, run, total_outflow, array_output, close_h5, get_const_info, get_model_info

Monte Carlo Runs
------------------

   .. method:: run_MC(nruns)

	When the model is constructed with the parameters according to the defined ranges, Monte Carlo simulations can
	be performed using the run_MC(number of runs) command

Re-Load an existing model structure
--------------------------------------

.. autofunction:: pyfuse.Load_model
