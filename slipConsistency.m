function [isConsist,dir_new] = slipConsistency(t,X,PHI,par,dir,dt)

dx=t*0; dd=dx;
for k=1:length(t)
    dx(k)=X(k,4);
    [~,dd(k)]=getd(t(k),X(k,:)',PHI,par);
end
dx2=dx+dd;

t_test=(t(end)-t(1))/10;
if t_test<dt
    t_test=0;
end
ind=(t>t(1)+t_test)&(t<(t(end)-t_test));
dx_short=dx(ind); dx2_short=dx2(ind);

% t_short=t(ind);
% bad_t1=t_short(dx_short~=dir(1)); bad_t2=t_short(dx2_short~=dir(1));
% (max(bad_t1)-min(bad_t1)<1e-3)&&(max(bad_t2)-min(bad_t2)<1e-3)

if all(sign(dx_short)==dir(1))&&all(sign(dx2_short)==dir(2))
    isConsist=1;
    dir_new=dir;
elseif (all(abs(dx)<1e-1)&&abs(dx(end))<1e-8)||(all(abs(dx2)<1e-1)&&abs(dx2(end))<1e-8)
    isConsist=-1;
    dir_new=dir;
else
    isConsist=0;
    dir_new=round([mean(sign(dx_short)) mean(sign(dx2_short))]);
end