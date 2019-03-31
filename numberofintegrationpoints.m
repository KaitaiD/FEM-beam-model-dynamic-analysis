function n = numberofintegrationpoints(ncoord,nelnodes,elident)

%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type

   if (ncoord == 1) 
     n = nelnodes;   
   elseif (ncoord == 2) 
     if (nelnodes == 3)
         n = 1;
     end
     if (nelnodes == 6)
         n = 3;
     end
     if (nelnodes == 4)
         n = 4;
     end
     if (nelnodes == 8)
         n = 9;
     end
   elseif (ncoord == 3) 
     if (nelnodes == 4)
         n = 1 ;
     end
     if (nelnodes == 10)
         n = 4;
     end
     if (nelnodes == 8)
         n = 8;
     end
     if (nelnodes == 20)
         n = 27;
     end
   end
end  
