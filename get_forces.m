function [ft1,fn1,ft2,fn2,tau] = get_forces(t,X,PHI,par,ie,dir)

mu_d1=par.mu_d1; mu_d2=par.mu_d2;

ft1=t*0; fn1=ft1; ft2=ft1; fn2=ft1; tau=zeros(2,length(t));
for k=1:length(t)
    tk=t(k);
    Xk=X(:,k);
    
    PHI_t=PHI(tk);
%     j1=PHI_t(1,1); j2=PHI_t(2,1);
    dj1=PHI_t(1,2); dj2=PHI_t(2,2);
    ddPhi=PHI_t(:,3);
    
    [Mb,B,G,W,dW,Eb] = Mat_inv(tk,Xk,PHI,par);
    
    if ie==1
        W1=W(1:2,:);
        dW1=dW(1:2,:);
        W2=W(3:4,:);
        wn2=W(4,:);
        dwn2=dW(4,:);
        
        L=[-mu_d2*dir;1];

        tmp=[Mb -W1' -W2'*L; [W1(:,1:3) ;wn2(:,1:3)] zeros(3,5)]\...
            ([-Eb;-W1(:,4:5);-wn2(:,4:5)]*ddPhi+[-B-G; -[dW1;dwn2]*[Xk(4:6);dj1;dj2]]);

        ft1(k)=tmp(6); fn1(k)=tmp(7);
        fn2(k)=tmp(8); ft2(k)=L(1)*fn2(k);
        
        if nargout>4
            tau(:,k)=tmp(4:5);
        end
        
    elseif ie==2
        wn1=W(2,:);
        dwn1=dW(2,:);
        W1=W(1:2,:);
        W2=W(3:4,:);
        dW2=dW(3:4,:);
        L=[-mu_d1*dir;1];

        tmp=[Mb -W1'*L -W2'; [wn1(:,1:3); W2(:,1:3)] zeros(3,5)]\...
            ([-Eb;-wn1(:,4:5);-W2(:,4:5)]*ddPhi+[-B-G; -[dwn1;dW2]*[Xk(4:6);dj1;dj2]]);

        fn1(k)=tmp(6); ft1(k)=L(1)*fn1(k);
        ft2(k)=tmp(7); fn2(k)=tmp(8);
        
        if nargout>4
            tau(:,k)=tmp(4:5);
        end

    elseif ie==3
        W1=W(1:2,:);
        wn1=W(2,:);
        dwn1=dW(2,:);
        W2=W(3:4,:);
        wn2=W(4,:);
        dwn2=dW(4,:);
        L1=[-mu_d1*dir(1);1];
        L2=[-mu_d2*dir(2);1];

        tmp=[Mb -W1'*L1 -W2'*L2; wn1(:,1:3) zeros(1,4); wn2(:,1:3) zeros(1,4)]\...
            ([-Eb;-wn1(:,4:5);-wn2(:,4:5)]*ddPhi+[-B-G; -[dwn1;dwn2]*[Xk(4:6);dj1;dj2]]);

        fn1(k)=tmp(6); ft1(k)=L1(1)*fn1(k);
        fn2(k)=tmp(7);  ft2(k)=L2(1)*fn2(k);
        
        if nargout>4
            tau(:,k)=tmp(4:5);
        end
        
    end

end

