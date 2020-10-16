clear variables
clc

syms x1 y1 q1 j1 j2 dx1 dy1 dq1 dj1 dj2 real
syms l beta m0 gama1 gama2 J1 J2 J0 k g al real

%% State vector
q=[x1 y1 q1 j1 j2].';
dq=[dx1 dy1 dq1 dj1 dj2].';

%% Position
r=[x1-(1/2)*l*cos(j1-q1), y1+(1/2)*l*sin(j1-q1);...
x1-l*cos(j1-q1)+(1/2)*beta*l*cos(q1), y1+l*sin(j1-q1)+(1/2)*beta*l*sin(q1);...
x1-l*cos(j1-q1)+beta*l*cos(q1)-(1/2)*l*cos(j2+q1), y1+l*sin(j1-q1)+beta*l*sin(q1)-(1/2)*l*sin(j2+q1)];

%% Velocity
v=r*0;
for n=1:5
    v=v+diff(r,q(n))*dq(n);
end

%% Energies
mm=[gama1*m0 m0 gama2*m0];
T=0; U=0;
for n=1:3
    T=T+(1/2)*mm(n)*(v(n,1)^2+v(n,2)^2);
    U=U+mm(n)*g*(cos(al)*r(n,2)-sin(al)*r(n,1));
end
T=T+(1/2)*J1*(dj1-dq1)^2+(1/2)*J0*dq1^2+(1/2)*J2*(dj2+dq1)^2;
U=U+(1/2)*k*((q1+j1)^2+(j2-j1)^2);

%% Assumptions

%J1=m*l^2/12; J2=beta*m*(beta*l)^2/12; J3=m*l^2/12;
% k1=k; k2=k;
T=eval(eval(T)); U=eval(eval(U));

%% Matricies

M=sym(zeros(5)); B=sym(zeros(5,1)); G=B;
for i=1:5
    for j=1:5
        M(i,j)=simplify(diff(T,dq(i),dq(j)));
        B(i)=B(i)+diff(T,dq(i),q(j))*dq(j);
    end
    B(i)=B(i)-diff(T,q(i));
    B(i)=simplify(B(i));
    G(i)=simplify(diff(U,q(i)));
end

%% Forces
constr=[x1;...
    y1;...
    x1-l*cos(j1-q1)+beta*l*cos(q1)-l*cos(j2+q1);...
    y1+l*sin(j1-q1)+beta*l*sin(q1)-l*sin(j2+q1)];
W=sym(zeros(4,5));
for i=1:length(constr)
    for j=1:5
        W(i,j)=simplify(diff(constr(i),q(j)));
    end
end
W=eval(W);

dW=0*W;
for n=1:5
    dW=dW+diff(W,q(n))*dq(n);
end 

% Lam=simplify((W*(M\W.'))\(W*(M\(B+G))-dW*dq));
% Lam=eval(eval(Lam));

%%
E=[0 0;0 0;0 0; 1 0; 0 1];

%% 
% NOTE: Need to add k
matlabFunction(M,B,G,W,dW,E,'File','Mat','Vars',{q1,j1,j2,dq1,dj1,dj2,[l,beta,m0,gama1,gama2,J1,J2,J0,k,g,al]})
% matlabFunction(Lam,'File','forcesStick')

%% d, dd, ddd
d=-l*cos(j1-q1)+beta*l*cos(q1)-l*cos(j2+q1);
dd=d*0; ddd=dd;
for k=3:5
    dd=dd+simplify(diff(d,q(k)))*dq(k);
end
dd=simplify(dd);

syms ddq1 ddj1 ddj2
ddq=[0;0;ddq1;ddj1;ddj2];
for k=3:5
    ddd=ddd+simplify(diff(dd,q(k)))*dq(k)+simplify(diff(dd,dq(k)))*ddq(k);
end
ddd=simplify(ddd);

matlabFunction(d,dd,ddd,'File','getd','Vars',{q1,dq1,ddq1,j1,dj1,j2,dj2,[ddj1;ddj2],[l,beta,m,gama,g]})

%% q, dq, ddq
y2=l*sin(j1-q1)+beta*l*sin(q1)-l*sin(j2+q1);
dy2=sym(0);
for k=1:5
    dy2=dy2+simplify(diff(y2,q(k)))*dq(k);
end
dq_exp=solve(dy2==0,dq1);

syms ddx ddq1 ddj1 ddj2
ddq=[ddx;0;ddq1;ddj1;ddj2];
ddq_exp=sym(0);
for k=3:5
    ddq_exp=ddq_exp+simplify(diff(dq_exp,q(k)))*dq(k)+simplify(diff(dq_exp,dq(k)))*ddq(k);
end
ddq_exp=simplify(ddq_exp)