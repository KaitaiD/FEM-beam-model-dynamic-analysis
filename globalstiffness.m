function Stif = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs)
%
%====================== Assemble the global stiffness matrix =================

   Stif = zeros(ndof*nnode,ndof*nnode);
   lmncoord = zeros(ncoord,maxnodes);
   lmndof = zeros(ndof,maxnodes);

%   Loop over all the elements
   for lmn = 1:nelem

%   Extract coords of nodes, DOF for the current element
      for a = 1:nelnodes(lmn)
          % extract nodal coordinates of current element
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        % extract the degrees of freedom
        for i = 1:ndof
          lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
      end
    n = nelnodes(lmn);
    ident = elident(lmn);
    kel = elstif(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof);

%   Add the current element stiffness to the global stiffness
    for a = 1:nelnodes(lmn)
      for i = 1:ndof
        for b = 1:nelnodes(lmn)
          for k = 1:ndof
            rw = ndof*(connect(a,lmn)-1)+i;
            cl = ndof*(connect(b,lmn)-1)+k;
% Global stiffness matrix = previous value + current element stiffness
            Stif(rw,cl) = Stif(rw,cl) + kel(ndof*(a-1)+i,ndof*(b-1)+k);
          end
        end
      end
    end
   end
end