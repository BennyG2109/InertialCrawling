function [ie_new,dir_new] = dchange(te,PHI,par,dt,ie_old)
global err choice_mode
% The test is performed some little time after the event, in case \dot{d}=0
t_test=te+dt;

%If dt=0, can take Xe instead of X_test
[q11,dq11]=completeX(PHI,t_test,par);
% x, dx doesnt affect d
q1=[0;0;q11]; dq1=[0;0;dq11];
X_test=[q1; dq1];

[~,~,ddd] = getd(t_test,X_test,PHI,par);

flag1=testmode(X_test,t_test,PHI,par,1,sign(ddd));
flag2=testmode(X_test,t_test,PHI,par,2,-sign(ddd));

ie1=flag1(1); ie2=flag2(1);
% If inconsistent - flag=[0 dir]
% dir is the potential future sliding direction
if length(flag1)>1
    dir1=flag1(2);
end
if length(flag2)>1
    dir2=flag2(2);
end

%TESTMODE:
% 1 - Mode is consistent,
%-1 - Contact break,
% 0 - Inconsistent mode
if ((ie1==-1)&&(ie2~=1))||((ie2==-1)&&(ie1~=1))  %Contact break consistent
    ie_new=4; dir_new=0;
elseif (ie1==1)&&(ie2~=1) %Stick-slip
    ie_new=1; dir_new=sign(ddd);
elseif (ie1~=1)&&(ie2==1) %Slip-stick
    ie_new=2; dir_new=-sign(ddd);
elseif (ie1==0)&&(ie2==0)
        ie_new=3;
        dir_new=[dir1 dir2];
elseif (ie1==1)&&(ie2==1)
    if dt<10
        [ie_new,dir_new] = dchange(te,PHI,par,5*dt,ie_old);
    else
        switch choice_mode
            case 1
                ie_new=ie_old;
            case 2
                ie_new=(ie_old==1)*2+(ie_old==2)*1;
            case 3
                ie_new=randi(2);
            case 0
                ie_new=NaN();
        end
        if ie_new==1
            dir_new=sign(ddd);
        else
            dir_new=-sign(ddd);
        end
    end
else
    ie_new=NaN(); dir_new=NaN();
end


end