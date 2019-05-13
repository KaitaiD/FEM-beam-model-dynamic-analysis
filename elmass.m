

function mel = elmass(ncoord,ndof,nelnodes,elident,lumpedmass,coord,materialprops)

%
%================= ELEMENT MASS MATRIX ====================================
%
% Calculate Jacobian + it's determinate.  Then assemble lumped mass matrix.
%

%
%  Assemble the (consistent) element mass matrix
%
%    Arguments:
%
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      elident            Element identifier (not used here - for future enhancements!)
%      coords[i,a]        ith coord of ath node
%      materialprops      Material properties passed on to constitutive procedures
%      displacement[i,a]  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi[i,inpt]         ith local coord of integration point no. intpt
%      w[intpt]           weight for integration point no. intpt
%      N[a]               Shape function associated with ath node on element
%      dNdxi[a,i]         Derivative of ath shape function wrt ith local coord
%      dNdx[a,i]          Derivative of ath shape function wrt ith global coord
%      dxdxi[i,j]         Derivative of ith global coord wrt jth local coord
%      dxidx[i,j]         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian

   npoints = numberofintegrationpointsmass(ncoord,nelnodes,elident);
   xi = zeros(ncoord,1);
   dxdxi = zeros(ncoord,ncoord);
   mel = zeros(ndof*nelnodes,ndof*nelnodes);
%
%  Set up integration points and weights    
%
   xilist = integrationpoints(ncoord,nelnodes,npoints,elident);
   w = integrationweights(ncoord,nelnodes,npoints,elident);

%  Loop over the integration points
   for intpt = 1 : npoints

%     Compute shape functions and derivatives wrt local coords
      for i = 1 : ncoord
        xi(i) = xilist(i,intpt);
      end    
      N = shapefunctions(nelnodes,ncoord,elident,xi);
      dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi);
      
%     Compute the jacobian matrix. Jacobian is a measure of volume change 
%     produced by deformation
      for i = 1 : ncoord
        for j = 1 : ncoord
          dxdxi(i,j) = 0.;
          for a = 1 : nelnodes
              % Jacobian = previous value + global coord * shape function
            dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(a,j);
          end
        end
      end
      
%   Calculate determinate of Jacobian
      dt = det(dxdxi);

      % Extract mass (rho) from input file
      rho = massdensity(materialprops);

      % Calculate the element masses in a 'non-lumped' matrix. (non-diag)
      for a = 1 : nelnodes
        for b = 1 : nelnodes
           for i = 1 : ndof
              row = ndof*(a-1)+i;
              col = ndof*(b-1)+i;
% Element mass = previous value + mass * shape_functions * integration weights * jacobian_determinate
              mel(col,row) = mel(col,row) + ...
                                  rho*N(b)*N(a)*w(intpt)*dt;
           end
        end
      end
   end

% Evaluate a lumped mass matrix using the row sum method (diagional matrix)
  if (lumpedmass) 
    for a = 1 : nelnodes*ndof
       for b = 1 : nelnodes*ndof
          if (a ~= b)
            mel(a,a) = mel(a,a) + mel(a,b);
            mel(a,b) = 0.; % elements on non-diagonal part are all zeros
          end
       end
    end
  end
end
