
function M = globalmass(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,lumpedmass,nelnodes,connect,materialprops)
%
%====================== Assemble the global mass matrix =================
%
%   Assemble the global mass matrix

   M = zeros(ndof*nnode,ndof*nnode);
   lmncoord = zeros(ncoord,maxnodes);

%   Loop over all the elements
   for lmn =  1 : nelem

%   Extract coords of nodes, DOF for the current element
      for a =  1 : nelnodes(lmn)
        for i = 1 : ncoord
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
      end
    n = nelnodes(lmn);
    ident = elident(lmn);
    % Calculate the lumped mass matrix
    mel = elmass(ncoord,ndof,n,ident,lumpedmass,lmncoord,materialprops);

%   Add the current element mass to the global mass
    for a = 1 : nelnodes(lmn)
        for i = 1: ndof
        for b = 1 : nelnodes(lmn)
          for k = 1 : ndof
            rw = ndof*(connect(a,lmn)-1)+i;
            cl = ndof*(connect(b,lmn)-1)+k;
          % Global mass = previous value + current element mass
            M(rw,cl) = M(rw,cl) + mel(ndof*(a-1)+i,ndof*(b-1)+k);
          end
        end
        end
    end
   end
  
  end