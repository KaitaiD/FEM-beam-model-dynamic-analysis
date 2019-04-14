# FEM-beam-model-dynamic-analysis
2D/3D Finite Element (FEM) beam model using Newmark Method for dynamic analysis

## Aim
This code provides an example of 2D and 3D dynamic linear elastic FEM code using __Newmark Method__. It contains a beam mesh with force exerted at the designated point and reveals the dynamic responses on the beam. It shows how to solve the global system equations in the FEM models.

## Input file
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

__dloads(_i_, _j_)__:         List of element tractions

              dloads(1,j) Element number
                            
              dloads(2,j) face number
                            
              dloads(3,j), dloads(4,j), dloads(5,j) Components of traction
                           
## Run the program
In order to run the program, it is necessary to have an input profile as explained in the last section. Then change the fopen command below to point to the file. Also change the fopen command in the post-processing step (near the bottom of the program) to point to a suitable output file. In the end go to __FEM_2Dor3D_linelast_dynamic_newmark.m__ as it is the main script that simulates the response. The main script has all been set up and the size of the model and the exerted force shall all be defined in the input profile.

## Output 
The script is able to plot the following things:

1. The initial mesh of the model before going through deformed process. It can be used as a check of the mesh of the FE model.

![alt text](https://github.com/KaitaiD/FEM-beam-model-dynamic-analysis/blob/master/f1.jpg)

2. The animation of the movement of the beam.

![alt text](https://github.com/KaitaiD/FEM-beam-model-dynamic-analysis/blob/master/f2.jpg)

3. The end displacement of the beam as a function of time

![alt text](https://github.com/KaitaiD/FEM-beam-model-dynamic-analysis/blob/master/f3.jpg)

In conclusion, the model is able to demonstrate the dynamic analysis of FE mesh given a certain force. Users feel free to change the size of the mesh and the position of the force and observe the change in the result.
