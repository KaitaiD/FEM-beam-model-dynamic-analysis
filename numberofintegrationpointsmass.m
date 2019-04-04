
function n = numberofintegrationpointsmass(ncoord,nelnodes,elident)
 %
%====================== No. integration points  for mass matrices =====================
%
 
   if (ncoord == 1) 
     n = nelnodes;   
   elseif (ncoord == 2) 
     if (nelnodes == 3)
         n = 4;
     end
     if (nelnodes == 6)
         n = 7;
     end
     if (nelnodes == 4)
         n = 4;
     end
     if (nelnodes == 8)
         n = 9;
     end
   elseif (ncoord == 3) 
     if (nelnodes == 4)
         n = 4 ;
     end
     if (nelnodes == 10)
         n = 5;
     end
     if (nelnodes == 8)
         n = 27;
     end
     if (nelnodes == 20)
         n = 27;
     end
   end
end   