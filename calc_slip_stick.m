function dX=calc_slip_stick(t,X,PHI,par,dir)
mu_d1=par.mu_d1;
% q1=X(3); dq1=X(6);

PHI_t=PHI(t);
% j1=PHI_t(1,1); j2=PHI_t(2,1);
dj1=PHI_t(1,2); dj2=PHI_t(2,2);
ddPhi=PHI_t(:,3);

[Mb,B,G,W,dW,Eb] = Mat_inv(t,X,PHI,par);

wn1=W(2,:);
dwn1=dW(2,:);
W1=W(1:2,:);
W2=W(3:4,:);
dW2=dW(3:4,:);
L=[-mu_d1*dir;1];

tmp=[Mb -W1'*L -W2'; [wn1(:,1:3); W2(:,1:3)] zeros(3,5)]\...
    ([-Eb;-wn1(:,4:5);-W2(:,4:5)]*ddPhi+[-B-G; -[dwn1;dW2]*[X(4:6);dj1;dj2]]);

dX=[X(4:6);tmp(1:3)];
end