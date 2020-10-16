function [E,E_pos]=energy(t,X,PHI,par,ie,dir)

[~,~,~,~,tau] = get_forces(t,X,PHI,par,ie,dir);

for k=1:length(t)
    tk=t(k);
    Xk=X(:,k);
    
    PHI_t=PHI(tk);
    dj1(k)=PHI_t(1,2); dj2(k)=PHI_t(2,2);
end

P1=tau(1,:).*dj1; P2=tau(2,:).*dj2;
E=trapz(t,P1)+trapz(t,P2);

E_pos=trapz(t,P1.*(P1>=0))+trapz(t,P2.*(P2>=0));