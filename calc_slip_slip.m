function dX=calc_slip_slip(t,X,PHI,par,dir)
mu_d1=par.mu_d1; mu_d2=par.mu_d2;

% q1=X(3); dq1=X(6);

PHI_t=PHI(t);
% j1=PHI_t(1,1); j2=PHI_t(2,1);
dj1=PHI_t(1,2); dj2=PHI_t(2,2);
ddPhi=PHI_t(:,3);

[Mb,B,G,W,dW,Eb] = Mat_inv(t,X,PHI,par);
W1=W(1:2,:);
wn1=W(2,:);
dwn1=dW(2,:);
W2=W(3:4,:);
wn2=W(4,:);
dwn2=dW(4,:);
L1=[-mu_d1*dir(1);1];
L2=[-mu_d2*dir(2);1];

tmp=[Mb -W1'*L1 -W2'*L2; wn1(:,1:3) zeros(1,4); wn2(:,1:3) zeros(1,4)]\...
    ([-Eb;-wn1(:,4:5);-wn2(:,4:5)]*ddPhi+[-B-G; -[dwn1;dwn2]*[X(4:6);dj1;dj2]]);

dX=[X(4:6);tmp(1:3)];
end



