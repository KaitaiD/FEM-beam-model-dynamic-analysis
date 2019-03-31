%
%================= ELEMENT STIFFNESS MATRIX ================================
%
function kel = elstif(ncoord,ndof,nelnodes,elident,coord,materialprops,displacement)
%
%  Assemble the element stiffness
%
%    Arguments;
%
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      elident            Element identifier (not used here - for future enhancements!)
%      coords(i,a)        ith coord of ath node
%      materialprops      Material properties passed on:constitutive procedures
%      displacement(i,a)  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi(i,inpt)         ith local coord of integration point no. intpt
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      strain(i,j)        strain_ij components
%      dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
%      kel(row,col)       Rows && cols of element stiffness
%
%
   npoints = numberofintegrationpoints(ncoord,nelnodes,elident);
   dNdx = zeros(nelnodes,ncoord);
   dxdxi = zeros(ncoord,ncoord);
   strain = zeros(ndof,ncoord);
   kel = zeros(ndof*nelnodes,ndof*nelnodes);
%
%  Set up integration points && weights    
%
   xilist = integrationpoints(ncoord,nelnodes,npoints,elident);
   w = integrationweights(ncoord,nelnodes,npoints,elident);
%
%  Loop over the integration points
%
   for intpt = 1:npoints

%     Compute shape functions && derivatives wrt local coords
%
      for i = 1:ncoord
        xi(i) = xilist(i,intpt);
      end      
      N = shapefunctions(nelnodes,ncoord,elident,xi);
      dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi);

      
%
%     Compute the jacobian matrix && its determinant
%
      for i = 1:ncoord
        for j = 1:ncoord
          dxdxi(i,j) = 0.;
          for a = 1:nelnodes
            dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(a,j);
          end
        end
      end
      
      dxidx = inv(dxdxi);
      dt = det(dxdxi);
%
%     Convert shape function derivatives:derivatives wrt global coords
%
      for a = 1:nelnodes
        for i = 1:ncoord
          dNdx(a,i) = 0.;
          for j = 1:ncoord
            dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
          end
        end
      end
%
%     Compute the (infinitesimal) strain by differentiating displacements
%     This step is not really necessary for linear elasticity calculations
%     where stiffness is independent of strain.  It is included:allow
%     extension:nonlinear materials later.
%
      for i = 1:ncoord
         for j = 1:ncoord
            strain(i,j) = 0.;
            for a = 1:nelnodes
              strain(i,j) = strain(i,j) + 0.5*(displacement(i,a)*dNdx(a,j)+displacement(j,a)*dNdx(a,i));
            end
         end
      end
%
%     Compute the material tangent stiffness (d stress/d strain)
%     ds/de is just C_ijkl for linear elasticity - this notation is used
%     to allow extension to nonlinear problems
%
      dsde = materialstiffness(ndof,ncoord,strain,materialprops);
%
%     Compute the element stiffness
%             
      for a = 1:nelnodes
        for i = 1:ndof
          for b = 1:nelnodes
            for k = 1:ndof
              row = ndof*(a-1)+i;
              col = ndof*(b-1)+k;
              for j = 1:ncoord
                for l = 1:ncoord
                  kel(col,row) = kel(col,row) + dsde(i,j,k,l)*dNdx(b,l)*dNdx(a,j)*w(intpt)*dt;
                end
              end
            end
          end
        end
      end
   end

end
