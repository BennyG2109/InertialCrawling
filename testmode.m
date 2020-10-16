function flag = testmode(X,t,PHI,par,ie,dir)
global err
% ie=1 - Stick-Slip
% ie=2 - Slip-Stick

%TESTMODE:
% 1 - Mode is consistent,
%-1 - Contact break,
% 0 - Inconsistent mode

% If inconsistent - flag=[0 dir]
% dir is the potential future sliding direction

mu_s1=par.mu_s1; mu_s2=par.mu_s2;

switch ie
    case 1      %Stick-Slip  
        [ft_stick,fn_stick,~,fn_slip] = get_forces(t,X,PHI,par,1,dir);
        mu_s=mu_s1;
    case 2      %Slip-Stick 
        [~,fn_slip,ft_stick,fn_stick] = get_forces(t,X,PHI,par,2,dir);
        mu_s=mu_s2;
end

if (fn_slip<=0)||(fn_stick<=0)   %Contact break
    flag=-1;
elseif (abs(ft_stick)/(mu_s*fn_stick)>=1+err)   %Inconsistent
    if ft_stick/(mu_s*fn_stick)>1+err
        flag=[0 -1];
    else
        flag=[0 1];
    end
else    %Consistent
    flag=1;
end

