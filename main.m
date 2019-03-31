%
%          Example 2D and 3D Linear elastic FEM code
%            Currently coded to run either plane strain or plane stress (2DOF) or general 3D but
%            could easily be modified for axisymmetry too.
%
%        Variables read from input file;
%        nprops              No. material parameters
%        materialprops(i)    List of material parameters
%        ncoord              No. spatial coords (2 for 2D, 3 for 3D)
%        ndof                No. degrees of freedom per node (2 for 2D, 3 for 3D)
%                            (here ndof=ncoord, but the program allows them to be different
%                            to allow extension to plate & beam elements with C^1 continuity)
%        nnode               No. nodes
%        coords(i,j)         ith coord of jth node, for i=1..ncoord; j=1..nnode
%        nelem               No. elements
%        maxnodes            Max no. nodes on any one element (used for array dimensioning)
%        nelnodes(i)         No. nodes on the ith element
%        elident(i)          An integer identifier for the ith element.  Not used
%                            in this code but could be used to switch on reduced integration,
%                            etc.
%        connect(i,j)        List of nodes on the jth element
%        nfix                Total no. prescribed displacements
%        fixnodes(i,j)       List of prescribed displacements at nodes
%                            fixnodes(1,j) Node number
%                            fixnodes(2,j) Displacement component number (1, 2 or 3)
%                            fixnodes(3,j) Value of the displacement
%        ndload              Total no. element faces subjected to tractions
%        dloads(i,j)         List of element tractions
%                            dloads(1,j) Element number
%                            dloads(2,j) face number
%                            dloads(3,j), dloads(4,j), dloads(5,j) Components of traction
%                            (assumed uniform) 
%
%  To run the program you first need to set up an input file, as described in 
%  the lecture notes.  Then change the fopen command below to point to the file.
%  Also change the fopen command in the post-processing step (near the bottom of the
%  program) to point to a suitable output file.  Then execute the file in
%  the usual fashion (e.g. hit the green arrow at the top of the MATLAB
%  editor window)
% 
%
% ==================== Read data from the input file ===========================
%
%
% YOU NEED TO CHANGE THE PATH & FILE NAME TO POINT TO YOUR INPUT FILE
%

clear all;
clc;

%infile=fopen('linear_elastic_Brick20.txt','r');
infile=fopen('soil.txt','r');
outfile=fopen('FEM_results.txt','w');

[nprops,materialprops,ncoord,ndof,nnode,coords,nelem,maxnodes,connect,nelnodes,elident,nfix,fixnodes,ndload,dloads] = read_input_file(infile);

fclose(infile);

% Plot the initial mesh as a check
close all
figure
plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,'g');

%
%============================ MAIN FEM ANALYSIS PROCEDURE ========================
%
%   dofs        Nodal displacements.  Let u_i^a be ith displacement component
%               at jth node.  Then dofs contain (u_1^1, u_2^1, u_1^2, u_2^2....) for 2D
%               and (u_1^1, u_2^1, u_3^1, u_1^2, u_2^2, u_3^2....) for 3D
%
%   K           Global stiffness matrix.  Stored as (K_1111 K_1112  K_1121  K_1122...
%                                                    K_1211 K_1212  K_1221  K_1222...
%                                                    K_2111 K_2112  K_2121  K_2122...)
%               for 2D problem and similarly for 3D problem
%   r           Force vector.  Currently only includes contribution from tractions
%               acting on element faces (i.e. body forces are neglected)
%
  dofs = zeros(ndof*nnode,1);  

  K = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs);

  r = globaltraction(ncoord,ndof,nnode,ndload,coords,nelnodes,elident,connect,dloads,dofs);

%
%  Fix constrained nodes.  We should really do this in a way that preserves the symmetry of K
%  but it's not worth it for the small demo problems here.
%
  for n = 1:nfix
     rw = ndof*(fixnodes(1,n)-1) + fixnodes(2,n);
     for cl = 1:ndof*nnode
        K(rw,cl) = 0;
     end
     K(rw,rw) = 1.;
     r(rw) = fixnodes(3,n);
  end
%
% Solve for the displacements
%

  dofs = K\r;

%
%================================= POST-PROCESSING =================================
%
% Create a plot of the deformed mesh
%

  defcoords = zeros(ndof,nnode); 
  scalefactor = 1.0;
  for i = 1:nnode
    for j = 1:ndof
       defcoords(j,i) = coords(j,i) + scalefactor*dofs(ndof*(i-1)+j); 
    end
  end

 figure
 plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,'g');
 hold on
 plotmesh(defcoords,ncoord,nnode,connect,nelem,elident,nelnodes,'r');  
  
 print_results(outfile, ...
    nprops,materialprops,ncoord,ndof,nnode,coords, ...
    nelem,maxnodes,connect,nelnodes,elident, ...
    nfix,fixnodes,ndload,dloads,dofs)
 
fclose(outfile);