function set = construct(A,B,v)
[V,J]=jordan(A);
% V=real(V);
% J=real(J);
V=eye(3);
J=A;

[m,n]=size(J);
b=inv(eye(m)-abs(J))*abs(inv(V)*B)*v;
A_=inv(V);
A=[A_;-A_];
b=[b;b];
set=Polyhedron(A,b);
end

