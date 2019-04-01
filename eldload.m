
function r = eldload(ncoord,ndof,nfacenodes,elident,coords,traction)
%
%====================== ELEMENT DISTRIBUTED LOAD VECTOR ==============
%
% Calculate number of integration points
  npoints = numberofintegrationpoints(ncoord-1,nfacenodes);

% Initialise matrices
  xi = zeros(ncoord-1,1);
  dxdxi = zeros(ncoord,ncoord-1);
  r = zeros(ndof*nfacenodes,1);

% Calculate integration points and weights
  xilist = integrationpoints(ncoord-1,nfacenodes,npoints);
  w = integrationweights(ncoord-1,nfacenodes,npoints);

% Loop through integration points
  for intpt = 1:npoints

    for i = 1:ncoord-1
      xi(i) = xilist(i,intpt);
    end

% Calculate shape functions and shape function derivatives
    N = shapefunctions(nfacenodes,ncoord-1,elident,xi);
    dNdxi = shapefunctionderivs(nfacenodes,ncoord-1,elident,xi);

%   Compute the jacobian matrix
    for i = 1:ncoord
      for j = 1:ncoord-1
        dxdxi(i,j) = 0.;
        for a = 1:nfacenodes
          dxdxi(i,j) = dxdxi(i,j) + coords(i,a)*dNdxi(a,j);
        end
      end
    end
    
    % Calculate Jacobian determinate depending on 2D or 3D
    if (ncoord == 2) 
      dt = sqrt(dxdxi(1,1)^2+dxdxi(2,1)^2);
    elseif (ncoord == 3) 
      dt = sqrt( ((dxdxi(2,1)*dxdxi(3,2))-(dxdxi(2,2)*dxdxi(3,1)))^2 ...
          + ((dxdxi(1,1)*dxdxi(3,2))-(dxdxi(1,2)*dxdxi(3,1)))^2 ...
          + ((dxdxi(1,1)*dxdxi(2,2))-(dxdxi(1,2)*dxdxi(2,1)))^2 );
    end
    
    % Calculate the total force/traction
   for a = 1:nfacenodes
      for i = 1:ndof
        row = ndof*(a-1)+i;
        % Force = previous value + shape function * traction * integration weights * jocobian determinate
        r(row) = r(row) + N(a)*traction(i)*w(intpt)*dt;
      end
    end
  end
end


