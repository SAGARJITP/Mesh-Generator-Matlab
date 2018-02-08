 function t1=checkinterior1(V,p)            
            if side(p,V(3,:),V(1,:),V(2,:))==1 && side(p, V(1,:),V(2,:),V(3,:))==1 && side(p, V(2,:),V(3,:),V(1,:))==1
               t1=1;
            else
               t1=0;
            end           
 end
        
 function t2=side(p1,p2,a,b)
            cp1= cross(b-a,p1-a);
            cp2= cross(b-a,p2-a);            
            
            if dot(cp1,cp2)>=0
                  t2=1;
            else
                  t2=0;
            end            
 end
        