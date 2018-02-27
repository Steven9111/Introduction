clear all
close all
clc

% PREDATOR PREY MODEL
%% constants
r=0.48;
c=0.01;
d=0.24;
e=0.005;

Nini = 100;
Wini=25;

%% Numerical solution 
totT = 100;
dt=0.001;
time=dt:dt:totT;

N(1) = Nini;
W(1)=Wini;

for t=2:numel(time)
        dN = (r*N(t-1)-c*W(t-1)*N(t-1))*dt;
        N(t) = N(t-1)+dN;
        
        dW = (-d*W(t-1)+e*W(t-1)*N(t-1))*dt;
        W(t) = W(t-1)+dW;
end

H = d*log(Nini)-e*Nini+r*log(Wini)-c*Wini;
H = d*log(N)-e*N+r*log(W)-c*W;

%Figure
figure
subplot(2,1,1)
plot(time,N)
hold on
plot(time, W)
legend('N','W')
hold off
xlabel('Time'),ylabel('Population')

subplot(2,1,2)
plot(W,N)

%% 
dt = linspace(0.0001,0.01,10);
dt = flip(dt);

totT = 100;
N(1) = Nini;
W(1)=Wini;

for i=1:numel(dt)
    time=dt(i):dt(i):totT;
    
for t=2:numel(time)
        dN = (r*N(t-1)-c*W(t-1)*N(t-1))*dt(i);
        N(i,t) = N(t-1)+dN;
        
        dW = (-d*W(t-1)+e*W(t-1)*N(t-1))*dt(i);
        W(i,t) = W(t-1)+dW;
        
        H = d*log(N)-e*N+r*log(W)-c*W;
                                                                                                                                                     
        %RMSE(i) = RMSE(N,
end
end

%% RMSE 
% At what timestep is the RMSE between the analytical and the numerical
% solution smaller than 0.001% of the analytical solition of the combined
% function (H)
%RMSE;    



%% Solve ODE with matlab solver:
tspan = [0 100];
y0 = [Wini, Nini];

[t,W_N] = ode45(@(t,W_N) ode_test(t,W_N,d,e,r,c), tspan, y0);

figure
plot(t,W_N(:,1),t,W_N(:,2))





