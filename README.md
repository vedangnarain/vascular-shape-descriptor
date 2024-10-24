# vascular-shape-descriptor

This is code to investigate the potential of a new shape descriptor to predict hypoxia. 

The networks are generated randomly as Voronoi tesselations. Each tesselation, called a 'selection', is used to produce smaller networks. Specifically, vessels are deleted in order of length, with the 'kill' count representing the number of vessels deleted, until flow is disconnected.

Here is how you can access a network:

1. Clone the repo.
2. Place `data.tar.gz` (ask me for it if you don't already have it) in the folder.
3. Use `visualise_simulation.py` for an example of how to access and plot the networks.  
4. Ignore `convert_oxygen_units.py`, unless you want to change the oxygen units.
