function [S,flag,Eout,Xx,Tt]=simu(par,om,PHI,cy_num,pl,rp)
global err choice_mode
S=[]; flag=[];

%flag=1 - good step
%flag=0 - fails cause take off or hit ground
%flag=-1 - numerical fail (bad forces, inconsistent, etc.)

choice_mode=0;
%choice mode
%0 - return NaN
%1 - keep old mode
%2 - switch mode
%3 - random mode

%pl:
%0 - no plot
%1 - animation
%2 - graphs
%3 - comparison between different simu
%4 - video

%% Simulation parameters
T=2*pi/om;
tf=cy_num*T;  %ending time

if nargin<5
    pl=0;   %1 - animation, 2 - graphs
end
if nargin<6
    rp=0;   %report? 1 or 0
end

dt=1e-6; %2*pi/om/20;
%For low freq. err=8e-3, for high freq. err=1e-4
if om<1
    err=10e-4;
else
    err=10e-4;
end

if pl==2
    col=get(gca,'Colororder');
    %close all
end

%% Simulation
% Completes the state vector
[q10,dq10]=completeX(PHI,0,par);
q0=[0;0;q10]; dq0=[0;0;dq10];
X0=[q0; dq0];

% Determine initial mode
[ie,dir] = dchange(0,PHI,par,dt,1);
% Fix leg velocity if slip
if ie==2
    [~,dx0]=getd(0,X0,PHI,par);
    X0(4)=-dx0;
end
% If hit ground or contact lost at initial point or numeric fail
if (ie==4)||(ie==5)
    warning('Houston we have a take off!')
    S=NaN();
    flag=0;
elseif (ie==6)||(ie==7)
    warning('We hit bottom')
    S=NaN();
    flag=0;
elseif isnan(ie)
    warning('Cant determine mode')
    S=NaN();
    flag=-1;
end
%Testing
% figure(5)
% plotter(0,X0',PHI,par,skip,1)

te=0; t=te; Xe=X0; ie_old=ie;
Tt=[]; Xx=[]; ft1=[]; fn1=[]; ft2=[]; fn2=[]; tau=[];
Etot=0; Etot_pos=0;

while (t(end)~=tf)&&all(ie~=[4:7])&&(isempty(S))
    % ie=1 - Stick-Slip
    % ie=2 - Slip-Stick
    % ie=3 - Slip-slip
    % ie=4/5 - Contact break
    % ie=6/7 - Another point hits the ground
    switch ie

    %%%%%%%%Stick-Slip
    case 1              
        if rp   %REPORT
            disp(['Mode stick-slip, direction ' num2str(dir)])
        end

        te_old=te;
        options=odeset('Events',@(t,X) event_stick_slip(t,X,PHI,par,dir),'reltol',1e-10,'abstol',1e-10);
        [t,X,te,Xe,ie]=ode45(@(t,X) calc_stick_slip(t,X,PHI,par,dir),[te tf],X0,options);

        %Save forces
        [ft1_i,fn1_i,ft2_i,fn2_i,tau_i] = get_forces(t,X',PHI,par,1,dir);
        if pl==2        
            figure(3)
            plot(t,ft1_i./fn1_i,'b',t,ft2_i./fn2_i,'k')
            hold on
            figure(5)
            plot(t,tau_i(1,:),'b',t,tau_i(2,:),'k')
            hold on
        elseif pl==3
            figure(1)
            semilogx([om om],[t(1) t(end)]/T,'Color',[col(4,:) 0.7],'linewidth',6)
            hold on
        end
        %Energy consumption
        [E,E_pos]=energy(t,X',PHI,par,1,dir);

        %%%MODE CHOICE
        if t(end)~=tf
            %If more than one mode detected, trying to clean numerics
            [te,Xe,ie]=clean_event(te,Xe,ie,te_old,dt,rp);

            %\dot{d}=0
            if all(ie==2)   
                [ie,dir]=dchange(te,PHI,par,dt,1);
                if isnan(ie)
                    warning('Oh no... Cant choose how to continue after d_dot=0')
                    S=NaN();
                    flag=-1;
                    break
                end
            %Transition to slip-slip
            elseif  all(ie==3)
                dir=[-1 dir];
            elseif all(ie==8)
                ie=3;
                dir=[1 dir];
            end
        end
        if ie==3, ie_old=1; end    %Can't figure why I did this

        if rp   %REPORT
            disp(['End time ' num2str(te') ' event ' num2str(ie')])
        end

    %%%%%%%%Slip-Stick
    case 2
        if rp   %REPORT
            disp(['Mode slip-stick, direction ' num2str(dir)])
        end

        te_old=te;
        options=odeset('Events',@(t,X) event_slip_stick(t,X,PHI,par,dir),'reltol',1e-10,'abstol',1e-10);
        [t,X,te,Xe,ie]=ode45(@(t,X) calc_slip_stick(t,X,PHI,par,dir),[te tf],X0,options);

        %Save forces   
        [ft1_i,fn1_i,ft2_i,fn2_i,tau_i] = get_forces(t,X',PHI,par,2,dir);
        if pl==2        
            figure(3)
            plot(t,ft1_i./fn1_i,'b',t,ft2_i./fn2_i,'k')
            hold on
            figure(5)
            plot(t,tau_i(1,:),'b',t,tau_i(2,:),'k')
            hold on
        elseif pl==3
            figure(1)
            semilogx([om om],[t(1) t(end)]/T,'Color',[col(1,:) 0.7],'linewidth',6)
            hold on
        end
        %Energy consumption
        [E,E_pos]=energy(t,X',PHI,par,2,dir);

        %%%MODE CHOICE   
        if t(end)~=tf
            %If more than one mode detected, trying to clean numerics
            [te,Xe,ie]=clean_event(te,Xe,ie,te_old,dt,rp);

            %\dot{x}=0
            if all(ie==1)
                [ie,dir]=dchange(te,PHI,par,dt,2);
                if isnan(ie)
                    warning('Oh no... Cant choose how to continue after d_dot=0')
                    S=NaN();
                    flag=-1;
                    break
                end
            %Transition to slip-slip
            elseif all(ie==3)
                dir=[dir -1];
            elseif all(ie==8)
                ie=3;
                dir=[dir 1];
            end
        end
        if ie==3, ie_old=2; end    %Can't figure why I did this

        if rp   %REPORT
            disp(['End time ' num2str(te') ' event ' num2str(ie')])
        end

    %%%%%%%%Slip-Slip  
    case 3
        if rp   %REPORT
            disp(['Mode slip-slip, directions ' num2str(dir)])
        end

        te_old=te; Xe_old=Xe;
        options=odeset('Events',@(t,X) event_slip_slip(t,X,PHI,par,dir),'reltol',1e-10,'abstol',1e-10);
        [t,X,te,Xe,ie]=ode45(@(t,X) calc_slip_slip(t,X,PHI,par,dir),[te tf],X0,options);

        isConsist=slipConsistency(t,X,PHI,par,dir,dt);
        if isConsist==0
            if err<0.05
                te=te_old; Xe=Xe_old; t=te;
                err=err*2;
                [ie,dir] = dchange(te,PHI,par,dt,ie_old);
                if isnan(ie)
                    warning('Cant determine mode')
                    S=NaN();
                    flag=-1;
                    break
                end
            else
                warning('Inconsistent slip-slip direction');       
                S=NaN();
                flag=-1;
                break
            end
        elseif isConsist==-1
            [te,Xe,ie]=clean_event(te,Xe,ie,te_old,dt,rp);
            ie=3;
        end

        %Save forces
        if (isConsist~=0)
            if om<1
                err=10e-4;
            else
                err=10e-4;
            end
            [ft1_i,fn1_i,ft2_i,fn2_i,tau_i] = get_forces(t,X',PHI,par,3,dir);
            if pl==2
                figure(3)
                plot(t,ft1_i./fn1_i,'b',t,ft2_i./fn2_i,'k')
                hold on
                figure(5)
                plot(t,tau_i(1,:),'b',t,tau_i(2,:),'k')
                hold on
            elseif pl==3
                figure(1)
                semilogx([om om],[t(1) t(end)]/T,'Color',[col(1,:) 0.7],'linewidth',6)
                hold on
                semilogx([om om],[t(1) t(end)]/T,'Color',[col(4,:) 0.5],'linewidth',6)
            end
            %Energy consumption
            [E,E_pos]=energy(t,X',PHI,par,3,dir);
        end
        
        %%%MODE CHOICE
        if (t(end)~=tf)&&(isConsist==1)
            %If more than one mode detected, trying to clean numerics
            [te,Xe,ie]=clean_event(te,Xe,ie,te_old,dt,rp);

            if all(ie==1)
                flag1=testmode(Xe',te,PHI,par,1,dir(2));
                if flag1(1)==1
                    dir=dir(2);
                else
                    ie=3;
                    dir(1)=-dir(1);
                end
            end
            if all(ie==2)
                flag2=testmode(Xe',te,PHI,par,2,dir(1));
                if flag2(1)==1
                    dir=dir(1);
                else
                    ie=3;
                    dir(2)=-dir(2);
                end
            end
        end

        if rp&&isConsist   %REPORT
            disp(['End time ' num2str(te') ' event ' num2str(ie')])
        end          
    end %end of switch

    
    %%%%Contact break or a link hitting the ground
    if (t(end)~=tf)
        if (ie==4)||(ie==5)
            warning('Houston we have a take off!')
            S=NaN();
            flag=0;
            break
        elseif (ie==6)||(ie==7)
            warning('We hit bottom')
            S=NaN();
            flag=0;
            break
        end
    else
        te=tf;
    end

    X0=Xe;
    if (abs(te-te_old)>1e-6)
        Xx=[Xx;X]; Tt=[Tt;t];
        ft1=[ft1;ft1_i]; fn1=[fn1;fn1_i]; ft2=[ft2;ft2_i]; fn2=[fn2;fn2_i];
        tau=[tau tau_i];
        Etot=Etot+E;
        Etot_pos=Etot_pos+E_pos;
    end

end

%%%Post processing
if ~isempty(Tt)
    plotter(Tt,Xx,PHI,par,pl,om,cy_num);
    if pl==2
        figure(3)
        xlabel('$t$ [sec]');
        ylabel('$f_{t,i}/f_{n,i}$');
        xlim([0 max(Tt)]); ylim([-1.1 1.1]);
        legend('$x_1$','$x_2$')
        hold on
        plot(Tt,Tt*0+par.mu_s1,'--',Tt,Tt*0-par.mu_s1,'--','Color',[0.6 0.6 0.6],'linewidth',2)
        plot(Tt,Tt*0+par.mu_d1,'--',Tt,Tt*0-par.mu_d1,'--','Color',[0.6 0.6 0.6],'linewidth',2)
        plot(Tt,Tt*0+par.mu_s2,'--',Tt,Tt*0-par.mu_s2,'--','Color',[0.6 0.6 0.6],'linewidth',2)
        plot(Tt,Tt*0+par.mu_d2,'--',Tt,Tt*0-par.mu_d2,'--','Color',[0.6 0.6 0.6],'linewidth',2)
        legend off
        figure(5)
        xlabel('$t$ [sec]');
        ylabel('Torques $\tau_i$');
        legend('$\tau_1$','$\tau_2$')
        xlim([0 max(Tt)]);
    elseif pl==3
        figure(1)
        xlabel('$\omega$ [rad/sec]');
        ylabel('$t/T$');
        ylim([2 3])
        set(gca,'ytick',[2 2.25 2.5 2.75 3],'yticklabel',{'$0$','$T/4$','$T/2$','$3T/4$','$T$'},'TickLabelInterpreter','latex')
    end
end

%% Analyze a single step
ind1=find(abs(ft1./fn1)>par.mu_s1+0.05); ind2=find(abs(ft2./fn2)>par.mu_s2+0.05);
if length(ind1)>2*cy_num||length(ind2)>2*cy_num
    warning('Bad forces');
    S=NaN();
    flag=-1;
end

if isempty(S)
    if choice_mode==2
        S=Xx(end,1)/(tf/T);
    else
        ind_halfstep=find(Tt>=tf-pi/om,1,'First');
        k=1; Sk=zeros(length(Tt)-ind_halfstep+1,1);
        for ind=length(Tt):-1:ind_halfstep
            ind_previous=find(Tt>=Tt(ind)-T,1,'First');
            Sk(k)=Xx(ind,1)-Xx(ind_previous,1);
            k=k+1;
        end
        S=mean(Sk); flag=1;
        if mad(Sk)>1e-3
            warning(['Step length might be wrong for omega=',num2str(om)])
        end
    end

    if nargout>2
        Eout=Etot_pos;
    end

    %if nargout>3
        %ind_last=find(Tt>tf-2*pi/om);
        %[~,dd] = getd(Tt(ind_last),Xx(ind_last,:)',PHI,par);
        %dx1=Xx(ind_last,4); dx2=dx1+dd;
        %slipslip=dx1((abs(dx1)>1e-6)&(abs(dx2)>1e-6));
        %slipslip_p=length(slipslip)/length(ind_last);
    %end
else
    if nargout>2
        Eout=NaN();
    end

end



end