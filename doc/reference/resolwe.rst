.. py:currentmodule:: orangecontrib.bio.resolwe

============================================================
Resolwe (:mod:`GenAPI`)
============================================================

Example of GenAPI usage:

.. literalinclude:: code/exampleGenesis.py

::

	Experiment id: 564a509e6b13390ffb40d4c8 - D. purpureum
	Experiment id: 564a54af6b13398f1640d4cd - Filter development
	Experiment id: 564a586b6b13398f1640d4d4 - cAMP pulses
	Experiment id: 56a944016b13395571175e36 - D. fasciculatum (Illumina)
	Experiment id: 56a944936b133909d8175e61 - D. lacteum (Illumina)
	Experiment id: 56a945fc6b133903da175e47 - D. fasciculatum (454)
	Experiment id: 56a946976b133909d8175e62 - P. pallidum (454)
	Experiment id: 56e2f39b6b1339964c33e7f0 - KO - gtaC knockout
	Experiment id: 56e2f50a6b1339964c33e7f1 - Wild Type
	Experiment id: 56e2f5a66b1339964c33e7f2 - CS - Cysteine-substituted strain
	Experiment id: 56e2f6706b1339964c33e7f3 - CM - Complemented mutant
	Experiment id: 56e2fdcc6b1339353d33e7e6 - D. discoideum

Download data from server:

.. literalinclude:: code/exampleGenesis1.py

::

	Time Points [0, 4, 8, 12, 16, 20, 24]
	Gene DPU_G0054708 [6.92022056529, 4.234260142395, 1.23427395022, 1.086050235081, 3.34525253324, 4.469964743405001, 3.26381321097]

Convert Json format to Orange data.Table:

.. literalinclude:: code/exampleGenesis2.py

::

	[71.582, 57.228, 57.609, 69.851, 65.566, 47.757, 27.165] {DPU_G0064618},
	 [181.716, 118.011, 100.540, 115.761, 53.397, 64.811, 56.353] {DPU_G0075396},
	 [3.382, 8.086, 1.581, 0.239, 2.893, 1.192, 1.781] {DPU_G0062146},
	 [0.758, 1.304, 0.678, 0.438, 1.509, 0.812, 0.751] {DPU_G0069396},
	 [629.349, 1290.830, 1493.227, 457.539, 15.799, 33.559, 114.367] {DPU_G0061882},
	 [0.665, 3.038, 2.136, 3.498, 5.774, 5.249, 5.878] {DPU_G0053074},
	 [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000] {DPU_G0055934},
	 [13.860, 43.443, 121.232, 272.804, 254.022, 156.572, 339.583] {DPU_G0053488},
	 [0.553, 0.495, 10.662, 35.259, 14.985, 22.432, 51.782] {DPU_G0070984},
	 [68.098, 55.382, 13.394, 13.612, 23.304, 17.440, 20.316] {DPU_G0062718}

.. autoclass:: orangecontrib.bio.resolwe.GenAPI
    :members: __init__, fetch_etc_objects, download_etc_data, etc_to_table
    :member-order: bysource
