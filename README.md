# FEM-beam-model-dynamic-analysis
2D/3D Finite Element (FEM) beam model using Newmark Method for dynamic analysis

### Aim
This code provides an example of 2D and 3D dynamic linear elastic FEM code using __Newmark Method__. It contains a beam mesh with force exerted at the designated point and reveals the dynamic responses on the beam. It shows how to solve the global system equations in the FEM models.

### Input file
In order to run the model, it is necessary to have the following variables as input:

__nprops__: No. material parameters

__materialprops(i)__: List of material parameters

__ncoord__:  No. spatial coords (2 for 2D, 3 for 3D)

__ndof__: No. degrees of freedom per node (2 for 2D, 3 for 3D) (here ndof=ncoord, but the program allows them to be different to allow extension to plate & beam elements with C^1 continuity)

__nnode__: No. nodes

__coords(i, j)__: ith coord of jth node, for i=1..ncoord; j=1..nnode

__nelem__:  No. elements
