# Cell_Evolution_stickymoves
source code for the results of the article: 
Evolution of multicellularity by collective integration of spatial information. 
For publication see the eLife: https://elifesciences.org/articles/56349

Based on the cellular potts model library TST, inherits requirement for Qt libraries (tested with Qt3 and Qt4) 

*** THIS BRANCH IS THE MOST UP TO DATE, and used to generate results for the publication ***

The code in its current form runs a simulation of cell evolution 
to investigate the circumstances under which grouping (sticking) of cells occurs. 
Cell adhesion is governed by a key-lock mechanism that evolves, and the amount of sticking at any given moment can be regulated
by the same factors as growth and movement in an evolvable manner.
Cells move around by persistant random walk and chemotaxis. 

Cells can eat food that is trickled in at random positions, either homogeneously or in patches.
Cells need food to maintain their mass, grow and divide, and to move. Division of resources among these tasks is evolvable and 
may be regulated through cell contact, cell size and amount of resources available.

The software requires the Qt libraries to compile (version 4 or higher). It's been extensively tested on Linux Ubuntu and other Linux distributions. 
Once compiled run with: 
Usage is: 
./cell_evolution path/to/parameters.par [optional arguments]

For a list of optional arguments type (note that you must always provide a valid path to a parameter file):
./cell_evolution path/to/parameters.par ARGUMENTS
