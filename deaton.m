% Shai Knight-Winnig 2016
% Deaton paper replication
%

clear
global alpha delta beta sigma sigma_sd mu CARA A
 


sigma=2;
sigma_sd=10;
alpha = 1;
delta = 0;
rho = 0.1;
beta  = 1/(1+rho);
r = 0.05;
R = 1+r;
mu=100;
CARA  =  1;
  
% grid size and vector of ones and number of iterations for simulation.
gridn   = 1000;
e       = ones(gridn,1);
  
% create grid for wealth
klow    = 20; khigh   = 320;
k = (klow: (khigh-klow)/(gridn-1) : khigh)';
 
%create grid for assets
alow = 0; ahigh = 200;
a = (alow: (ahigh-alow)/(gridn-1) : ahigh)';
 
%"Benchmark" consumption
css=60;
 
% consumption matrix (today's k in column, tomorrow's in rows)
C = e*k' - a*e';

% make negative consumption unattractive
if sigma ~= 'exp'
    C = max(C , 0);         % take out negative consumption (or else u(c) would not be defined)
    U = u(C) + ~C*(-1e50);  % add penalty for zero consumption
else
    U = u(C);
end
 
%Generate random numbers used for computing expectation
y=normrnd(mu,sigma_sd,gridn,1); % we generate lots of random numbers to take empirical expectations
X=R*a*e'+e*y'; % we construct the matrix (1+r)*a+y for all possible combinations of a and y.
 
 
%initialize iteration
crit    = 1; iter = 1; 
v  = u(k - k/R)/(1-beta);      % initial guess (value of keeping k and c constant) in the stochastic setting
%v =0*v + 5;                    % alternative initial guess
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate on Bellman equation  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while crit > 1e-5
    
    V=interp1(k,v(:,iter),X); % we calculate the value function at the matrix X
    EV=mean(V,2); % we take empirical expectations
    vnew = max( U + beta*EV*e' )'; % the Bellman step.
 
    % compute normalized criterion (scale is in consumption per period)
    if sigma ~= 'exp'
        critnew = (1-beta)*max(abs(vnew - v(:,iter)))/css^(-sigma+1);
    else
        critnew = (1-beta)*max(abs(vnew - v(:,iter)))/css/exp(-CARA*css);
    end
    
    % display iteration information
    % note: contraction mapping implies that critnew/crit < beta
    display([iter , critnew, critnew/crit ])
    crit = critnew; iter = iter +1; v(: , iter)  = vnew;
    
end
% compute policy function

V=interp1(k,v(:,iter),X);
EV=mean(V,2);
[vnew, policykindex]= max( U + beta*EV*e');
policyk = a(policykindex);
figure(1);
hold on;
plot(k,k,k,k-policyk,k,policyk)
axis([0 200 0 200])
title('Optimal Consumption in Income Fluctuation Problem')
legend('45-Degree Line','Optimal Consumption','Optimal Asset Holding','Location','NorthWest')
set(gca,'YTick',0:20:140)
set(gca,'XTick',0:40:280)
xlabel('wealth')
ylabel('consumption')  
hold off;

% simulate
a=zeros(1000,500);
income=zeros(1000,500);


for j=1:1000
    y=normrnd(mu,sigma_sd,500,1);

wealth=zeros(1,500)';
asset(1)=20;
wealth(1)=(1+r)*asset(1)+y(1);
consumption=zeros(1,500)';
asset=zeros(1,500)';
time=linspace(1,500,500)';

for i=1:500
    wealth(i)=(1+r)*asset(i)+y(i);
    if i<500
        asset(i+1)=policyk(floor((wealth(i)-20)/0.3)+1);
        consumption(i)=wealth(i)-asset(i+1);
    end
end
consumption(500)=wealth(500);  
a(j,:)=asset';
income(j,:)=y';
end

figure(2);
subplot(2,3,1)
hist(a(:,2));
title('t=2')
subplot(2,3,2)
hist(a(:,10));
title('t=10')
subplot(2,3,3)
hist(a(:,50));
title('t=50')
subplot(2,3,4)
hist(a(:,100));
title('t=100')
subplot(2,3,5)
hist(a(:,500));
title('t=500')

abar=mean(a,1);
incomebar=mean(income,1);
ratio=abar./incomebar;
figure(3);
plot(time,ratio);
title('ratio of average asset over average income')





