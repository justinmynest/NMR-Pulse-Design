function [Mdot]= rfevolve(t,M,dw,A,V,n)
Ox=[0 0 0; 0 0 -1;0 1 0];
Oy=[0 0 1; 0 0 0;-1 0 0];
Oz=[0 -1 0; 1 0 0;0 0 0];

u = 1; w = 0;

if n>2
    i = n; j = n;
    while j>1
        
        u = 2*cos(V(j)*t)*u;
        j = j-1;
    end
    u = A(i)*u;
    i = i-1;    
    while i>1        
        v = 1; k = i-1;        
        while k>1            
            v = v*2*cos(V(k)*t);
            k = k-1;
        end
        v = v*2*A(i)*sin(V(i)*t);
        w = w + v;
        i = i-1;
    end
    w = w + u;
else
    w = 2*A(2)*cos(V(2)*t);                     % earlier done with A(1) !! check
end


Mdot=(dw*Oz+A(1)*Ox+w*Oy)*M;
end


