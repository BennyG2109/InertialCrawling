function [eventvalue,stopthecalc,eventdir] = event_stick_slip(t,X,PHI,par,dir)
l1=par.l; l2=par.l*par.beta;
mu_s1=par.mu_s1;

q1=X(3);
PHI_t=PHI(t);
j1=PHI_t(1,1);

[ft1,fn1,~,fn2] = get_forces(t,X,PHI,par,1,dir);

% [W1;wn2]*[X(6:8);dj1;dj2]

[~,dd,~] = getd(t,X,PHI,par);
hit=[l1*sin(pi-(j1-q1)) l1*sin(pi-(j1-q1))+l2*sin(q1)];

% ie=1 - Stick-Slip
% ie=2 - Slip-Stick
% ie=3 - Slip-slip
% ie=4/5 - Contact break
% ie=6/7 - Another point hits the ground
eventvalue=[1,dd,ft1-mu_s1*fn1,fn1,fn2,hit,ft1+mu_s1*fn1]';
stopthecalc=[0,1,1,1,1,1,1,1]';
eventdir=[0,0,1,-1,-1,-1,-1,-1]';

end