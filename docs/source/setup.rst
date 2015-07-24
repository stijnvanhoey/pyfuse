Model setup functions
==================================

Model structure construction options
-------------------------------------   

.. autofunction:: pyFUSE.Set_options

Template file for model options::

	######################################################################
	##    Model options input file
	##    The options defined here are used to set up the model structure
	##
	##    (1) upper-layer architecture: uplayer 
	##        tension1_1: upper layer broken up into tension and free storage
	##        tension2_1: tension storage sub-divided into recharge and excess
	##        onestate_1: upper layer defined by a single state variable
	##        surface1_1: upper layer defined by a surface storage representing and a tension reservoir
	##    (2) lower-layer architecture and baseflow: lowlayer_baseflow
	##        tens2pll_2: tension reservoir plus two parallel tanks
	##        unlimfrc_2: baseflow resvr of unlimited size, frac rate
	##        unlimpow_2: baseflow resvr of unlimited size, power recession
	##        fixedsiz_2: baseflow reservoir of fixed size
	##    (3) surface runoff: surface
	##        arno_x_vic: ARNO/Xzang/VIC parameterization (upper zone control)
	##        prms_varnt: PRMS variant (fraction of upper tension storage)
	##        tmdl_param: TOPMODEL parameterization (only valid for TOPMODEL qb)
	##        oflwtresh: threshold based overland flow generation
	##    (4) percolation     
	##        perc_f2sat: water from (field cap to sat) avail for percolation
	##        perc_w2sat: water from (wilt pt to sat) avail for percolation
	##        perc_lower: perc defined by moisture content in lower layer (SAC)
	##        perc_tresh: threshold based percolation
	##        perc_nodrain: percolation represents the baseflow routing
	##    (5) evaporation
	##        sequential: sequential evaporation model
	##        rootweight: root weighting
	##    (6) interflow   
	##        intflwnone: no interflow
	##        intflwsome: linear interflow
	##        intflwtresh: threshold based interflow generation   
	##    (7) routing
	##        rout_all1: touting combined subflows
	##        no_rout: no routing
	##        rout_ind: rout subflows individual
	######################################################################
	## MODEL DECISION 	=   OPTION
	uplayer				=	onestate_1
	lowlayer_baseflow	= 	unlimfrc_2
	surface				=	arno_x_vic 
	percolation			=	perc_w2sat
	evaporation			=	sequential  
	interflow			=	intflwnone
	routing				=	rout_all1

Some well-known model structures can be selected by calling the name of the model,

for the NAM model structure [1]_::

	>>> pf.Set_options(filename=False,default='NAM')
	NAM model options are selected
	The selected options are
	{'evaporation': 'sequential',
	 'interflow': 'intflwtresh',
	 'lowlayer_baseflow': 'unlimfrc_2',
	 'percolation': 'perc_tresh',
	 'reservoirs': '220',
	 'routing': 'rout_ind',
	 'surface': 'oflwtresh',
	 'uplayer': 'surface1_1'}

for the PDM model structure [2]_::

	>>> pf.Set_options(filename=False,default='PDM')
	PDM model options are selected
	The selected options are
	{'evaporation': 'sequential',
	 'interflow': 'intflwnone',
	 'lowlayer_baseflow': 'unlimfrc_2',
	 'percolation': 'perc_w2sat',
	 'routing': 'rout_ind',
	 'reservoirs': '200',
	 'surface': 'arno_x_vic',
	 'uplayer': 'onestate_1'}


	
Model parameters
------------------

.. autofunction:: pyFUSE.Set_pars

.. autofunction:: pyFUSE.Set_pars_for_run	

Template file for parameters::

	###############################################
	##    Model Parameter input file
	##    The parameter is defined by his distribution, boundaries and extra info needed by distribution
	##    provide on each line one parameter with following information:
	##
	##    name : string
	##        Name of the parameter
	##    minval : float
	##        Minimum value of the parameter distribution
	##    maxval :  float
	##        Maximum value of the parameter distribution
	##    optguess : float
	##        Optimal guess of the parameter, must be between min and max value
	##    pardistribution : string
	##        choose a distributionfrom: randomUniform, randomTriangular, randomTrapezoidal, randomNormal, randomLogNormal
	##    *kargs  : de
	##        Extra arguments necessary for the chosen distribution
	#################################################
	## NAME	MIN MAX OPTGUESS DISTRIBUTION ARGS*
	S1max	50. 5000.000 400. randomTriangular 1000.
	S2max	100. 10000.000 1000. randomNormal 500. 25. 
	fitens 0.01 1.0 0.99 randomLogNormal	0.5 0.2
	firchr 0.050 0.950 0.5 randomTrapezoidal 0.4 0.6
	fibase 0.050 0.950 0.5 randomUniform
	r1 0.050 0.950 0.5 randomUniform
	ku 0.01 1000. 0.044 randomUniform
	c 0.99 20.0 1. randomUniform
	alfa 1.000 250. 150. randomUniform
	psi 1.000 5.0 2.5 randomUniform
	kappa 0.050 0.950 0.5 randomUniform
	ki 0.001 1000. 0.00833 randomUniform
	ks 0.001 10000. 0.5 randomUniform
	n 1.000 10. 3. randomUniform
	v 0.00001 0.250 0.004	randomUniform
	vA 0.001 0.250	0.0015 randomUniform
	vB 0.001 0.250 0.0015 randomUniform
	Acmax 0.050 0.950 0.5 randomUniform
	b 0.001 3.0 0.2 randomUniform
	loglambda	5.000 10.0 7.5 randomUniform
	chi 2.000 5.0 3.5 randomUniform
	mut 0.010 5.0 0.6 randomUniform
	be 0.99 4. 3.1 randomUniform
	alfah	0.01 0.99 0.5 randomUniform
	tg 0.0 0.7 0.3 randomUniform
	tif 0.0 0.7 0.26 randomUniform
	tof 0.0 0.7 0.12 randomUniform
	ko 0.01 0.99 0.15 randomUniform
	timeo 2. 48 24 randomUniform
	timei 2. 250. 20 randomUniform
	timeb 200. 10000. 2100. randomUniform

	
Model constant values
----------------------
	
.. autofunction:: pyFUSE.Set_cnts

Template file for constants::

	#######################################################################
	##	Model constants values
	##	The options defined here are used to set up the model structure
	##		
	##	  area: the size of the catchment in km2
	##	  timestep: timestep relative to the hourly timestep 	
	#######################################################################	
	## CONSTANT			 =   VALUE
	area				=	362.
	timestep			= 	1.


	
References
^^^^^^^^^^^^

.. [1] DHI. MIKE 11, A Modelling System for Rivers and Channels, Reference Manual. Horsholm, Denmark: DHI Water & Environment, 2008.
.. [2] Moore, R J. The PDM Rainfall-runoff Model. Hydrology and Earth System Sciences 11, no. 1 (2007): 483–499.