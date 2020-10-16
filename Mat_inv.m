function [Mb,B,G,W,dW,Eb] = Mat_inv(t,X,PHI,par)
q1=X(3); dq1=X(6);
PHI_t=PHI(t);
j1=PHI_t(1,1); j2=PHI_t(2,1); dj1=PHI_t(1,2); dj2=PHI_t(2,2);

[M,B,G,W,dW,E] = Mat(q1,j1,j2,dq1,dj1,dj2,par);

Mb=[M(:,1:3) -E];
Eb=M(:,4:5);

end

