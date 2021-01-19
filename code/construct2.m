function set = construct2(A,B,v,c)
center=inv(eye(3)-A)*B(:,10:12)*c;
[V,J]=jordan(A);
% V=real(V);
% J=real(J);

V=eye(3);
J=A;

[m,n]=size(J);
b=inv(eye(m)-abs(J))*(abs(inv(V)*B)*v);
A_=inv(V);
A=[A_;-A_];
b=[b;b];
set=Polyhedron(A,b)+center;
end

