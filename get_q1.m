function q1=get_q1(j1,j2,par)
l1=par.l; l2=par.l*par.beta;

%express r12 in frame attached to L2, then rotate by -q1
r1=-l2/2+l1*exp(-j1*1j);
r2=l2/2+l1*exp((pi+j2)*1j);
r21=r2-r1;
q1=-angle(r21);

end