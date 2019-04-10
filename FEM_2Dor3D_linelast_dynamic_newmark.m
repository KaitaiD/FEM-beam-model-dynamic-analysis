% DYNAMIC linear elastic code for elastodynamics

% This brings together all the other functions to calculate
% the overall solution in terms of displacement

clear all
clc
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
infile=fopen('Linear_elastic_dynamic_beam.txt','r');
[nprops,materialprops,ncoord,ndof,nnode,coords,nelem,maxnodes,connect,nelnodes,elident,nfix,fixnodes,ndload,dloads] = read_input_file(infile);
fclose(infile);

% Plot the initial mesh as a check (undeformed shape)
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
lumpedmass = true;

% Define the GLOBAL mass, stiffness and traction matrices
%--------------------------------------------------------
M = globalmass(ncoord,ndof,nnode,coords, ...
    nelem,maxnodes,elident,lumpedmass,nelnodes,connect,materialprops);
K = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs);
F = globaltraction(ncoord,ndof,nnode,ndload,coords,nelnodes,elident,connect,dloads,dofs);

%  Fix constrained nodes (boundary conditions).
%--------------------------------------------------------
for n = 1: nfix
    rw = ndof*(fixnodes(1,n)-1) + fixnodes(2,n);
    for cl = 1 : ndof*nnode
        K(rw,cl) = 0;
        K(cl,rw) = 0.;
        M(rw,cl) = 0.;
        M(cl,rw) = 0.;
    end
    M(rw,rw) = 1.;
    F(rw) = fixnodes(3,n);
end

% Define Newmark parameters (original  = 0.5, 0) unconditionally stable
beta1 = 0.6;
beta2 = 1.01*((2*beta1 + 1)^2)/16;
termination_time = 150;           % seconds
timestep = 0.01;
nprint = 100;
count = 0;

dt = timestep;
total_dt = termination_time/timestep;
nsteps = total_dt;

% Initialise matricies
un = zeros(nnode*ndof,1);
vn = zeros(nnode*ndof,1);
an1 = zeros(nnode*ndof,1);
% Store time and end displacement for plotting
vxy = zeros(2,nsteps);

%    Newmark time stepping scheme (problem is linear so all matrices are constant)
%    Initial accelerations are -M^-1(Ku + F)

% Assuming no damping - calculate K (MK - 3.126)
MK = M + 0.5*beta2*dt*dt*K;
an = M\(-K*un + F);

%--------------------------------------------------------
% Commence time integration - Loop through every time step
%--------------------------------------------------------

for n = 1:nsteps
    % Calculate F_residual (3.127)
    rn = F - K*(un +dt*vn + 0.5*(1.-beta2)*dt*dt*an );
    rn = F - K*(un +timestep*vn + 0.5*timestep*timestep*an ); 
    
    an1 = (MK^-1)*rn;
    
    % Calculate the velocities and displacements based upon acceleration
    vn1 = (vn + dt*(1.-beta1)*an + dt*beta1*an1);
    un1 = (un + dt*vn + (1.-beta2)*0.5*dt*dt*an + 0.5*beta2*dt*dt*an1 );
    un = un1;
    vn = vn1;
    an = an1;
    vxy(1,n) = n*dt;
    vxy(2,n) = un(2*nnode);
    
    %      Plot displaced mesh every nprint steps.
    %--------------------------------------------------------
    if (count == nprint)
        scale_factor = 3.;
        defcoords = zeros(ndof,nnode);
        for i = 1 : nnode
            for j = 1 : ndof
                defcoords(j,i) = coords(j,i) + scale_factor*un(ndof*(i-1)+j);
            end
        end
        clf
        axis([0,11,-5.5,5.5])
        % Apply plotmesh function for animation
        plotmesh(defcoords,ncoord,nnode,connect,nelem,elident,nelnodes,'r');
        pause(0.025);
        count = 0;
    end
    count = count + 1;
    % End time step looping
end


% Plot the end displacement as a function of time (2D graph)
%--------------------------------------------------------
figure
plot(vxy(1,:),vxy(2,:))
xlabel({'Time'},'FontSize',16);
ylabel({'End displacement'},'FontSize',16);

