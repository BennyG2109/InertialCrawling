function [d,dd,ddd] = getd(t,X,PHI,par)
l1=par.l; l2=par.l*par.beta;

j1=0*t; j2=j1; dj1=j1; dj2=j1; ddj1=j1; ddj2=j1; q1=j1; dq1=j1;
for k=1:length(t)
    q1(k)=X(3,k); dq1(k)=X(6,k);
    
    tk=t(k);
    PHI_t=PHI(tk);
    j1(k)=PHI_t(1,1); j2(k)=PHI_t(2,1);
    dj1(k)=PHI_t(1,2); dj2(k)=PHI_t(2,2);
    
    if nargout>2
        ddj1(k)=PHI_t(1,3); ddj2(k)=PHI_t(2,3);
    end
end
    

d = -l1.*cos(j2+q1)-l1.*cos(j1-q1)+l2.*cos(q1);
if nargout > 1
    dd = -dq1.*(l1.*sin(j1-q1)-l1.*sin(j2+q1)+l2.*sin(q1))+dj1.*l1.*sin(j1-q1)+dj2.*l1.*sin(j2+q1);
end
if nargout > 2
    [~,~,ddq1]=completeX(PHI,t,par);
    ddd = dj2.^2.*l1.*cos(j2+q1)+dj1.^2.*l1.*cos(j1-q1)+ddj1.*l1.*sin(j1-q1)+ddj2.*l1.*sin(j2+q1)-ddq1.*l1.*sin(j1-q1)+ddq1.*l1.*sin(j2+q1)-ddq1.*l2.*sin(q1)+l1.*cos(j2+q1).*dq1.^2+l1.*dq1.^2.*cos(j1-q1)-l2.*dq1.^2.*cos(q1)+dj2.*dq1.*l1.*cos(j2+q1).*2.0-dj1.*dq1.*l1.*cos(j1-q1).*2.0;
end
