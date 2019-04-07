# FEM-beam-model-dynamic-analysis
2D/3D Finite Element (FEM) beam model using Newmark Method for dynamic analysis

### Aim
This code provides an example of 2D and 3D dynamic linear elastic FEM code using __Newmark Method__. It contains a beam mesh with force exerted at the designated point and reveals the dynamic responses on the beam. It shows how to solve the global system equations in the FEM models.

### Input file
In order to run the model, it is necessary to have the following variables as input:

__nprops__: No. material parameters

__materialprops(_i_)__: List of material parameters

__ncoord__:  No. spatial coords (2 for 2D, 3 for 3D)

__ndof__: No. degrees of freedom per node (2 for 2D, 3 for 3D) (here ndof=ncoord, but the program allows them to be different to allow extension to plate & beam elements with C^1 continuity)

__nnode__: No. nodes

__coords(_i_, _j_)__: _i_ th coord of _j_ th node, for _i_ =1..ncoord; _j_ =1..nnode

__nelem__:  No. elements

__maxnodes__:  Max no. nodes on any one element (used for array dimensioning)

__nelnodes(_i_)__:  No. nodes on the _i_ th element

__elident(_i_)__:  An integer identifier for the _i_ th element.  Not used in this code but could be used to switch on reduced integration, etc.

__connect(_i_, _j_)__:  List of nodes on the _j_ th element

__nfix__:    Total no. prescribed displacements

__fixnodes(_i_, _j_)__:       List of prescribed displacements at nodes

                         fixnodes(1,j) Node number
                            
                         fixnodes(2,j) Displacement component number (1, 2 or 3)
                            
                         fixnodes(3,j) Value of the displacement
                            
__ndload__:  Total no. element faces subjected to tractions

__dloads(i,j)__         List of element tractions

                        dloads(1,j) Element number
                            
                        dloads(2,j) face number
                            
                        dloads(3,j), dloads(4,j), dloads(5,j) Components of traction
                           
