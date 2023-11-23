function [Mdot]= evolveconstant(t,M,dw,Ax,Ay)
Ox=[0 0 0; 0 0 -1;0 1 0];
Oy=[0 0 1; 0 0 0;-1 0 0];
Oz=[0 -1 0; 1 0 0;0 0 0];
w=1;
v=1;
u=1;
Mdot=((dw*w)*Oz+(Ax*u)*Ox+(Ay*v)*Oy)*M;
end




