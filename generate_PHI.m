function PHI = generate_PHI(n,om,A1,A2,phi,j0,A1n,B1n,A2n,B2n)

syms ts

if length(j0)==1
    fun=[j0; j0];
else
    fun=[j0(1); j0(2)];
end

switch n
    %Regular sine, mostly for simulation
    case 1
        fun=fun+[A1*cos(om*ts+0.5*phi);A2*cos(om*ts-0.5*phi)];
        
    %Fourier series with constant phase    
    case 2
        A1n=[1 A1n]; B1n=[0 B1n]; %First harmony amp. is constant
        scale=sum(abs(A1n+B1n*1i));
        for l=1:length(A1n)
            fun=fun+[A1/scale*(A1n(l)*cos(l*(om*ts+0.5*phi))+B1n(l)*sin(l*(om*ts+0.5*phi)));...
                A2/scale*(A1n(l)*cos(l*(om*ts-0.5*phi))+B1n(l)*sin(l*(om*ts-0.5*phi)))];
        end
        
    %Fourier series of both angles (fitted phase)
    case 3
        scale1=sum(abs(A1n+B1n*1i)); scale2=sum(abs(A2n+B2n*1i));
        for l=1:length(A1n)
            fun=fun+[A1/scale1*(A1n(l)*cos(l*om*ts)+B1n(l)*sin(l*om*ts));...
                A2/scale2*(A2n(l)*cos(l*om*ts)+B2n(l)*sin(l*om*ts))];
        end
        
    %Time scaling (one angle)
    case 4
        ds=sym(0);
        for l=1:length(A1n)
            ds=ds+A1n(l)*cos(l*ts)+B1n(l)*sin(l*ts);
        end
        T=2*pi/om; tt=linspace(0,T,100);
        ds_t=eval(subs(ds,ts,tt));
        amp=min(ds_t);
        if amp<-1
            ds=ds/abs(amp);
        end
        ds=ds+1;
        
        s=int(ds,ts,[0,ts]);
        
        %My scaling:
        fun=fun+[A1*cos(subs(s,ts,om*ts+0.5*phi));A2*cos(subs(s,ts,om*ts-0.5*phi))];
        %Yizhar's scaling:
        %fun=fun+[A1*cos(subs(s,ts,om*ts)+0.5*phi);A2*cos(subs(s,ts,om*ts)-0.5*phi)];
        
    %Time scaling (both angles)
    case 5
        ds1=sym(1); ds2=sym(1);
        for l=1:length(A1n)
            ds1=ds1+A1n(l)*cos(l*om*ts)+B1n(l)*sin(l*om*ts);
            ds2=ds2+A2n(l)*cos(l*om*ts)+B2n(l)*sin(l*om*ts);
        end
        s1=int(om*ds1,ts,[0,ts]);
        s2=int(om*ds2,ts,[0,ts]);
        
        fun=fun+[A1*cos(s1);A2*cos(s2)];
end

PHI=matlabFunction([fun diff(fun) diff(fun,2)]);

end

