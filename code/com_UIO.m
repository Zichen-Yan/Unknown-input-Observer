function [H,T,K1,N,K2] = com_UIO(A,B,C,E1,M)
H = E1*inv((C*E1)'*(C*E1))*(C*E1)';
MT=B-H*C*B;
T=M\MT;
A1=A-H*C*A;
MK1 = (place(A1',C',[0.8,0.2,0.5]))';
K1=M\MK1;
N=M\((A-H*C*A-M*K1*C)*M);
K2=M\((A-H*C*A-M*K1*C)*H);
end

