% Decentralized State Feedback Controller
clc; clear; close all;

A=[1 1 0 1 0 1
  -1 0 -1 0 0 1
  1 0 0 -1 1 1
  -1 1 -1 0 0 0
  -1 -1 1 1 -1 -1
  0 -1 0 0 -1 0];
B1=[0 -1 -1
    0 0 0
    -1 -1 1
    -1 0 0
    0 0 1
    -1 1 1];
B2=[0  0  0
   -1  0  1
   -1  1  0
    1 -1  0
   -1  0 -1
    0  1  1];
C1=[0 1 0 -1 -1 -1
    0 0 0 -1  0  0
    1 0 0  0 -1  0];
D11=[0 0 1
    -1 0 0
     0 0 0];
D12=[0 1 1
    0 0 0
    1 1 1];
eta=0.00000001;

% LMIs
X=sdpvar(6,6);
Z=sdpvar(3,6,'full');
gamma=sdpvar(1);
Const=[];
Const=[Z(1,3:6)==0;Z(2:3,1:2)==0];
Const=[Const; X >= eta*eye(size(X)); X(1:2,3:6)==0; X(3:6,1:2)==0];
M=[X*A'+A*X+Z'*B2'+B2*Z         B1         X*C1'+Z'*D12'; 
    B1'                             -gamma*eye(3)                    D11';
    C1*X+D12*Z                   D11             -gamma*eye(3)];
Const=[Const; M <= eta*eye(size(M))];
optimize(Const,gamma);
H_infinity_gain=value(gamma)
X=value(X);
E_X=eig(X);
Z=value(Z);
K=Z*inv(X)
