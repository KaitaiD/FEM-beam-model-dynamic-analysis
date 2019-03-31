function C = materialstiffness(ndof,ncoord,strain,materialprops)

%================= Material Stiffness ==================================
%
%    Computes elasticity tensor C_{ijkl} = shear modulus and Poissons ratio
%    Currently coded either for plane strain, plane stress or general 3D.
%

   mu = materialprops(1);
   nu = materialprops(2);
   
   C = zeros(ndof,ncoord,ndof,ncoord);
   
   if (ncoord == 2)   

%  planestrain = 0 => plane stress, planestrain = 1 => plane strain
   planestrain = materialprops(3);
 
     for i = 1:2
       for j = 1:2
         for k = 1:2
           for l = 1:2
             if (planestrain==1) 
               if (i==j && k==l)
                   C(i,j,k,l) = C(i,j,k,l)+2*mu*nu/(1-2*nu);
               end
             else
               if (i==j && k==l)
                   C(i,j,k,l) = C(i,j,k,l)+2*mu*nu/(1-nu);
               end
             end
             if (i==l && k==j)
                 C(i,j,k,l) = C(i,j,k,l)+mu;
             end
             if (i==k && j==l)
                 C(i,j,k,l) = C(i,j,k,l)+mu;
             end
           end
         end
       end
     end
     
   elseif (ncoord == 3) 

     for i = 1:3 
       for j = 1:3
         for k = 1:3
           for l = 1:3
             if (i==j && k==l)
                 C(i,j,k,l)=C(i,j,k,l) + 2.*mu*nu/(1.-2.*nu);
             end
             if (i==k && j==l)
                 C(i,j,k,l)=C(i,j,k,l) + mu;
             end
             if (i==l && j==k)
                 C(i,j,k,l)=C(i,j,k,l) + mu;
             end
            end
          end
        end
      end
    end

 end
