function r = globaltraction(ncoord,ndof,nnodes,ndload,coords,nelnodes,elident,connect,dloads,dofs)

%
%===================== Assemble the global traction vector =============
%

   r = zeros(ndof*nnodes,1);
   traction = zeros(ndof,1);

   for load = 1:ndload
%
%     Extract the coords of the nodes on the appropriate element face
%
      lmn = dloads(1,load);
      face = dloads(2,load);
      n = nelnodes(lmn);
      ident = elident(lmn);
      nfnodes = nfacenodes(ncoord,n,ident,face); 
      nodelist = facenodes(ncoord,n,ident,face);     
      lmncoord = zeros(ncoord,nfnodes);
      for a = 1:nfnodes
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(nodelist(a),dloads(1,load)));
        end
        for i = 1:ndof
          lmndof(i,a) = dofs(ndof*(connect(nodelist(a),dloads(1,load))-1)+i);
        end
      end
%
%    Compute the element load vector
%
     for i = 1:ndof
       traction(i) = dloads(i+2,load);
     end

     rel = eldload(ncoord,ndof,nfnodes,ident,lmncoord,traction);
%
%    Assemble the element load vector into global vector
%
     for a = 1:nfnodes
       for i = 1:ndof
         rw = (connect(nodelist(a),dloads(1,load))-1)*ndof+i;
         r(rw) = r(rw) + rel((a-1)*ndof+i);
       end
     end

   end
end   

