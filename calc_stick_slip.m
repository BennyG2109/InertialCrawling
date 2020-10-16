function dX=calc_stick_slip(t,X,PHI,par,dir)
mu_d2=par.mu_d2;

% q1=X(3); dq1=X(6);

PHI_t=PHI(t);
% j1=PHI_t(1,1); j2=PHI_t(2,1);
dj1=PHI_t(1,2); dj2=PHI_t(2,2);
ddPhi=PHI_t(:,3);

[Mb,B,G,W,dW,Eb] = Mat_inv(t,X,PHI,par);

W1=W(1:2,:);
dW1=dW(1:2,:);
W2=W(3:4,:);
wn2=W(4,:);
dwn2=dW(4,:);
L=[-mu_d2*dir;1];

tmp=[Mb -W1' -W2'*L; [W1(:,1:3) ;wn2(:,1:3)] zeros(3,5)]\...
    ([-Eb;-W1(:,4:5);-wn2(:,4:5)]*ddPhi+[-B-G; -[dW1;dwn2]*[X(4:6);dj1;dj2]]);

dX=[X(4:6);tmp(1:3)];
end



