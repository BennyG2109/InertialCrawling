clc
close all
clear variables

%%
par.l=0.17; par.beta=1.1;
par.m0=0.194; par.gama1=0.1082; par.gama2=0.1082;
par.J1=9.87e-5; par.J2=9.87e-5; par.J0=7.763e-4;

par.al=deg2rad(0); %slope - not in use
par.g=9.81;

par.k=0; %stiffness - not in use

%Hard
mu_s=0.1723; mu_d=0.0895;
%Soft
%mu_s=0.3983; mu_d=0.0831;

par.mu_s1=mu_s; par.mu_d1=mu_s;
par.mu_s2=mu_s; par.mu_d2=mu_s;
%Can vary the friction coef. among legs
%and between dynamic _d and static _s coefs.

%% One run
A1=deg2rad(18); %Amplitude 1
A2=deg2rad(18); %Amplitude 2
j0=deg2rad(110); %Nominal angle \phi_0
phi=-deg2rad(20); %Phase
om=8; %Freq.
cy_num=4;

PHI = generate_PHI(1,om,A1,A2,phi,j0);
%[S,flag,Eout,Xx,Tt]=simu(par,om,PHI,cy_num,pl,rp)
[S1,flag,Eout,Xx,Tt]=simu(par,om,PHI,cy_num,2,0)

%flag=1 - good step
%flag=0 - fails cause take off or hit ground
%flag=-1 - numerical fail (bad forces, inconsistent, etc.)

%pl:
%0 - no plot
%1 - animation
%2 - graphs
%3 - comparison between different simu
%4 - video

%rp: 1 or 0 - plots reports of mode changes to cmd. window