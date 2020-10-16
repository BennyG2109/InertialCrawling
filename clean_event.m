function [te_new,Xe_new,ie_new] = clean_event(te,Xe,ie,te_old,dt,rp)

if length(te)>1
    while (abs(te(1)-te_old)<1e-3)&&(length(te)>1)
        te=te(2:end);
        ie=ie(2:end);
        Xe=Xe(2:end,:);
    end
end

if length(te)>1
    if te(1)-te(2)>dt
        if rp   %REPORT
            disp(['Chose time ' num2str(te(1)) ' event ' num2str(ie(1)) 'and not time ' num2str(te(2))])
        end
        ie=ie(1);
        te=te(1);
        Xe=Xe(1,:);
    end   
end
ie_new=ie;
te_new=te;
Xe_new=Xe;
end