function [eventvalue,stopthecalc,eventdir] = event_slip_stick(t,X,PHI,par,dir)
l1=par.l; l2=par.l*par.beta;
mu_s2=par.mu_s2;

q1=X(3);
PHI_t=PHI(t);
j1=PHI_t(1,1);

[~,fn1,ft2,fn2] = get_forces(t,X,PHI,par,2,dir);

dx=X(4);
% [~,dd,~] = getd(t,X,PHI,PAR);
hit=[l1*sin(pi-(j1-q1)) l1*sin(pi-(j1-q1))+l2*sin(q1)];

% ie=1 - Stick-Slip
% ie=2 - Slip-Stick
% ie=3 - Slip-slip
% ie=4/5 - Contact break
% ie=6/7 - Another point hits the ground
eventvalue=[dx,1,ft2-mu_s2*fn2,fn1,fn2,hit,ft2+mu_s2*fn2]';
stopthecalc=[1,0,1,1,1,1,1,1]';
eventdir=[0,0,1,-1,-1,-1,-1,-1]';
end