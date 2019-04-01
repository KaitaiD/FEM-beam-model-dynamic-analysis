
function list = facenodes(ncoord,nelnodes,elident,face)
%======================= Lists of nodes on element faces =============
%
%    This procedure returns the list of nodes on an element face
%    The nodes are ordered so that the element face forms either
%    a 1D line element or a 2D surface element for 2D or 3D problems

% This is a list of the node numbers from the input file

   i3 = [2,3,1];
   i4 = [2,3,4,1]; 

   list = zeros(nfacenodes(ncoord,nelnodes,face),1);

   if (ncoord == 2) 
     if (nelnodes == 3) 
       list(1) = face;
       list(2) = i3(face);
     elseif (nelnodes == 6) 
       list(1) = face;
       list(2) = i3(face);
       list(3) = face+3;
     elseif (nelnodes==4) 
       list(1) = face;
       list(2) = i4(face);
     elseif (nelnodes==8) 
       list(1) = face;
       list(2) = i4(face);
       list(3) = face+4;
     end
   elseif (ncoord == 3) 
     if (nelnodes==4) 
       if   (face == 1)
           list = [1,2,3];
       elseif (face == 2)
           list = [1,4,2];
       elseif (face == 3)
           list = [2,4,3];
       elseif (face == 4)
           list = [3,4,1];
       end
     elseif (nelnodes == 10) 
       if   (face == 1)
           list = [1,2,3,5,6,7];
       elseif (face == 2)
           list = [1,4,2,8,9,5];
       elseif (face == 3)
           list = [2,4,3,9,10,6];
       elseif (face == 4)
           list = [3,4,1,10,8,7];
       end
     elseif (nelnodes == 8) 
       if   (face == 1)
           list = [1,2,3,4];
       elseif (face == 2)
           list = [5,8,7,6];
       elseif (face == 3)
           list = [1,5,6,2];
       elseif (face == 4)
           list = [2,3,7,6];
       elseif (face == 5)
           list = [3,7,8,4];
       elseif (face == 6)
           list = [4,8,5,1];
       end
     elseif (nelnodes == 20)  
       if   (face == 1)
           list = [1,2,3,4,9,10,11,12];
       elseif (face == 2)
           list = [5,8,7,6,16,15,14,13];
       elseif (face == 3)
           list = [1,5,6,2,17,13,18,9];
       elseif (face == 4)
           list = [2,6,7,3,18,14,19,10];
       elseif (face == 5)
           list = [3,7,8,4,19,15,20,11];
       elseif (face == 6)
           list = [4,8,5,1,20,16,17,12];
       end
     end
   end
end