
function n = nfacenodes(ncoord,nelnodes,elident,face)
 %====================== No. nodes on element faces ================
%
%   This procedure returns the number of nodes on each element face
%   for various element types.  This info is needed for computing
%   the surface integrals associated with the element traction vector

if (ncoord == 2) 
% For linear 2D elements there are 2 nodes on each element face
     if (nelnodes == 3 || nelnodes == 4)
         n = 2;
% For quadratic 2D elements there are 3 nodes on each element face
     elseif (nelnodes == 6 || nelnodes == 8)
         n=3;
     end
   elseif (ncoord == 3) 
     if (nelnodes == 4)
         n = 3;
     elseif (nelnodes == 10)
         n = 6;
     elseif (nelnodes == 8)
         n = 4;
     elseif (nelnodes == 20)
         n = 8;
     end
   end
end