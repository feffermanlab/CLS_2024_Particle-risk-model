% This code runs a basic SIR model to track community spread
% Gives initial conditions for Multi_Room_Model.m
% So this code needs to be run before Multi_Room_Model.m
% To use this code for the IC, capture S, I, R values
    %i.e. run [t,S,I,R] = Community
    
%Graphical outputs: plots the S, I, R over time

%author: Anna Sisk
%email: annahsisk@gmail.com
%Last update: 05/11/2022
function [t,S,I,R] = Community(S0,I0,R0, MaxTime)
%
% 
% Sets up default parameters if necessary.
%inital conditions-currently arbitrary
if nargin == 0
  %S0=9000+randi(10,1,1);
 % I0=10000-S0;
 % R0=0;
  S0= 1;
  I0= 0.00000127;
  R0= 0;
  MaxTime=144; %length of outbreak
end
tvec = linspace(0, MaxTime, MaxTime+1);
% The main iteration 
options = odeset('RelTol', 1e-3);
[t, pop]=ode45(@SIR,tvec,[S0,I0,R0],options);
   S=pop(:,1); 
   I=pop(:,2);
   R=pop(:,3);

   figure; %plotting the solution curves
   %font size is 16
   plot(t,S,'-g',t,I,'-k',t,R,'-b','LineWidth',1.0);  
   legend('S','I','R')
   ylabel('Population')
   xlabel('Time')
   title('Long Term Population Behavior')
 
end
function dPop=SIR(t,pop, parameter)
%parameters-not currently based on data
 %beta=.00002;  %contract/transmission rate
 %gamma=.1; %recovery rate
beta = 6/10;
gamma = .33;

       
dPop=zeros(3,1);

S=pop(1); 
I=pop(2);
R=pop(3);

%model-equations (1) in the write up
dPop(1)=-beta*S*I;
dPop(2)=beta*S*I-gamma*I;
dPop(3)=gamma*I;
end