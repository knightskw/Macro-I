% Shai Knight-Winnig 2016

clear;

%initialize parameters:
beta=0.9;
alpha=0.3;
A=18.5;
N=100; %number of states/controls (can update here without making any other code changes)
time_periods=100; %number of time periods (can update here to visually test convergence without any other code changes needed)
ubd=20; %state/control space upper bound
lbd=ubd/N; %state/control space lower bound
dist=(ubd-lbd)/(N-1); %distance between state/control points
E=(log(A*(1-alpha*beta))+(beta*alpha/(1-alpha*beta))*log(A*beta*alpha))/(1-beta); %V*() theoretical solution
consF=alpha/(1-alpha*beta); %This is the constant F from theoretical solution

%NOTE: this problem assumes that you have the same number of states/controls
%and that you have a positive starting capital amount
%initialize vectors & matrices:
x=linspace(lbd,ubd,N); %evenly spaced states from (0,20], assuming initial capital > 0
y=linspace(lbd,ubd,N); %control (policy rules) - same as states in this problem
optValue=zeros(N,time_periods+1); %matrix with optimal values for each state/time period
optPolicy=zeros(N,time_periods);  %matrix with decision rules for optimal value at each state/time period

for t=time_periods:-1:1                             %iterate over time
    for i=1:length(x)                               %loop for every x (state)
        utility=zeros(N,1);                         %vector to track utility values
        for j=1:length(y)                           %calculate utility for every y (control)
            if y(j)<=A*(x(i)^alpha)                 %only evaluate y(j) if feasible
                index=cast(y(j)/dist,'uint8');      %gives index for control choice that we pass into next period
                utility(j) = log(A*(x(i))^alpha-y(j)) + beta*optValue(index,t+1);
            else
                utility(j) = -1000000;              %assign very low utility to infeasible controls    
            end
        end
    [M I] = max(utility);      %return max value (M) and index of max value (I) after iterating across all controls
    optValue(i,t)=M;           %store maximum in optValue matrix (for particular state/time period
    optPolicy(i,t)=y(I);       %store optimal policy rule for each state/time period
    clearvars utility;         %clear for iteration to next state
    end
end

%%PLOTS - UNCOMMENT AS NEEDED TO GENERATE%%

%plot simulated value function vs. v*() theoretical solution

plot(x,optValue(:,1),'b*'); %removing x=0 (state) from plot since utility is very low
hold on;
plot(x,E+consF*log(x),'r','linewidth',2)
xlabel('Initial State');
ylabel('Maximum Utility');
legend('V(x)','V*(x)');
title('Simulation V() versus V*()');
hold off;


%plot simulated decision rules vs. theoretical solution

plot(x,optPolicy(:,1),'b*');
hold on;
plot(x,alpha*beta*A*x.^alpha,'r','linewidth',2)
xlabel('Initial State');
ylabel('Optimal Decision (Capital to Save for Next State)');
legend('Simulated Decision Rule','Theoretical Decision Rule');
title('Simulated Policy versus Theoretical Optimal Policy');
hold off;
