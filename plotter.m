function []=plotter(t,X,PHI,par,pl,om,cy_num)
%pl:
%------------
%0 - no plot
%1 - animation
%2 - graphs
%3 - comparison between different simu
%4 - video
%------------

l1=par.l; l2=par.l*par.beta;

if nargin<7
    cy_num=1;
end

switch pl
    case 1 %1 - animation
        skip=3;
        dt=gradient(t/pi*om);
        tl=length(t);
        x_lim=[-4*l1+min(X(:,1)) max(X(:,1))+4*l1];
        for k=1:skip:tl
            x1=X(k,1);
            y1=X(k,2);
            q1=X(k,3);
            PHI_t=PHI(t(k));
            j1=PHI_t(1,1);
            j2=PHI_t(2,1);

            xx=[x1,x1+l1*cos(pi-(j1-q1)),x1+l1*cos(pi-(j1-q1))+l2*cos(q1),x1+l1*cos(pi-(j1-q1))+l2*cos(q1)+l1*cos(pi-(j2+q1))];
            yy=[y1,y1+l1*sin(pi-(j1-q1)),y1+l1*sin(pi-(j1-q1))+l2*sin(q1),y1+l1*sin(pi-(j1-q1))+l2*sin(q1)-l1*sin(pi-(j2+q1))];
            
            plot(xx,yy,'color',[0 0 k/tl])
            axis equal
            xlim(x_lim); ylim([-0.1*l1 3*l1]);
            hold on
            plot(x_lim,[0 0],'k','linewidth',1)
            
            % Slip direction arrow
            dx=X(k,4);
            dq1=X(k,6); dj1=PHI_t(1,2); dj2=PHI_t(2,2);
            [~,dd] = getd(t(k),X(k,:)',PHI,par);
            if abs(dx)>1e-5
                quiver(x1,-0.1*l1,sign(dx),0,0.01,'k')
            end
            if abs(dx+dd)>1e-5
                quiver(xx(end),-0.1*l1,sign(dd),0,0.01,'k')
            end

            shg
            pause(dt(k));
            hold off
        end
        
    case 2 %2 - graphs          
        x1=X(:,1); q1=X(:,3);
        j1=x1*0; j2=j1;
        for k=1:length(t)
            PHI_t=PHI(t(k));
            j1(k)=PHI_t(1,1); j2(k)=PHI_t(2,1);
        end
        
        m0= par.m0; m1=par.gama1*m0; m2=par.gama2*m0;
        [d,dd] = getd(t,X',PHI,par);
        %d=l1*cos(pi-(j1-q1))+l2*cos(q1)+l1*cos(pi-(j2+q1));
        x2=x1+d;
        xc=m1*(x1+0.5*l1*cos(pi-(j1-q1)))+m0*(x1+l1*cos(pi-(j1-q1))+0.5*l2*cos(q1))+m2*(x1+l1*cos(pi-(j1-q1))+l2*cos(q1)+0.5*l1*cos(pi-(j2+q1)));
        xc=xc/(m1+m2+m0);
        
        figure(1)
        %subplot(2,1,1)
        plot(t,j1*180/pi,t,j2*180/pi,'linewidth',2)
        y_lim=ylim;
        hold on
        plot((t(end)-2*pi/om)*[1 1],y_lim,'--k','linewidth',1.5)
        hold off
        xlim([0 max(t)]); ylim(y_lim);
        xlabel('$t$ [sec]'); ylabel('Angles [deg]')
        legend('$\varphi_1$','$\varphi_2$')
        grid on
        %subplot(3,1,2)
        %plot(j1,j2,'linewidth',2)
        %xlabel('$\varphi_1$'); ylabel('$\varphi_2$')
        %subplot(2,1,2)
        %plot(t,dd,'linewidth',2)
        %y_lim=ylim;
        %hold on
        %plot((t(end)-2*pi/om)*[1 1],y_lim,'--k','linewidth',1.5)
        %hold off
        %xlabel('$t$ [sec]'); ylabel('$\dot{d}$');
        %xlim([0 max(t)]); ylim(y_lim);

        xmid=x1+0.5*d;
        figure(2)
        %plot(t,xc-xmid,t,d,'linewidth',2)
        %hold on
        plot(t,x1,t,x2,'linewidth',2)
        y_lim=ylim;
        hold on
        plot((t(end)-2*pi/om)*[1 1],y_lim,'--k','linewidth',1.5)
        hold off
        %plot(t,t*0,'--k')
        xlim([0 max(t)]); ylim(y_lim);
        xlabel('t [sec]');
        %legend('$\Delta_c=x_c-x_{mid}$','$d(t)=x_2-x_1$','$x_1$','$x_2$')
        legend('$x_1$','$x_2$')
        grid on

        dx1=X(:,4);
        figure(4)
        plot(t,dx1,'b',t,dx1+dd,'k')
        y_lim=ylim;
        hold on
        plot((t(end)-2*pi/om)*[1 1],y_lim,'--k','linewidth',1.5)
        hold off
        xlim([0 max(t)]); ylim(y_lim);
        xlabel('$t$ [sec]'); ylabel('$\dot{x}_i$ [m/s]')
        legend('$\dot{x}_1$','$\dot{x}_2$')
        grid on
        
%     case 3
%         x1=X(:,1);
%         T=t/(2*pi/om);
%         figure(1)
%         plot(T,x1,'DisplayName',['$\omega=$',num2str(om)])
%         hold on
%         
% %         h=gca;
% %         p=h.Children;
% %         c=p(1).Color;
% %         i=round(length(t)/2);
% %         text(T(i),x1(i),['$\omega=$',num2str(om)],'EdgeColor',c,'Backgroundcolor','w','linewidth',2);
%         xlim([0 max(T)])
%         xlabel('$T=t\omega/2\pi$')
%         ylabel('$x_1$');
%         legend show

    case 4 %video
        hz=30*om/2/pi;
        v = VideoWriter('anim','MPEG-4');
        v.FrameRate=hz;
        open(v)
            
        x_lim=[-0.1*l1+min(X(:,1)) max(X(:,1))+3*l1];
        tot_frames=floor(cy_num*(2*pi/om)*hz);
        for t_k=1:tot_frames
            k=find(t>=t_k/hz,1,'First');
            
            x1=X(k,1);
            y1=X(k,2);
            q1=X(k,3);
            PHI_t=PHI(t(k));
            j1=PHI_t(1,1);
            j2=PHI_t(2,1);

            xx=[x1,x1+l1*cos(pi-(j1-q1)),x1+l1*cos(pi-(j1-q1))+l2*cos(q1),x1+l1*cos(pi-(j1-q1))+l2*cos(q1)+l1*cos(pi-(j2+q1))];
            yy=[y1,y1+l1*sin(pi-(j1-q1)),y1+l1*sin(pi-(j1-q1))+l2*sin(q1),y1+l1*sin(pi-(j1-q1))+l2*sin(q1)-l1*sin(pi-(j2+q1))];
            
            plot(xx,yy,'color',[t_k/tot_frames 0 1])
            axis equal
            y_lim=[-0.1*l1 2*l1];
            xlim(x_lim); ylim(y_lim);
            hold on
            plot(x_lim,[0 0],'k','linewidth',1)
            plot([0 0],y_lim,':','color',0.4*[1 1 1]);
            
            % Slip direction arrow
            dx=X(k,4);
            dq1=X(k,6); dj1=PHI_t(1,2); dj2=PHI_t(2,2);
            [~,dd] = getd(t(k),X(k,:)',PHI,par);
            if abs(dx)>1e-5
                quiver(x1,-0.05*l1,sign(dx),0,0.05,'k')
            end
            if abs(dx+dd)>1e-5
                quiver(xx(end),-0.05*l1,sign(dd),0,0.05,'k')
            end
            set(gca,'XTick',[],'YTick',[])
            shg
            pause(1/hz);
            hold off
            
            frame = getframe(gca);
            writeVideo(v,frame);
        end
        
        close(v);
end