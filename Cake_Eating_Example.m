% Shai Knight-Winnig 2016
% Cake eating problem,T=2 


clear;
%Define Grid and parameters 
x=linspace(.1,10,100);
y=linspace(0.1,10,100);
v1=zeros(size(x));
v2=zeros(size(x));
v3=zeros(size(x));
A=18.5
alpha=.3
beta=.9;

%Given v3, solve for v2 
for i=1:length(i)                               %loop for every x
    for j=1:i                                   %calculate rhs for every y
        rhs(j) = log(x(i)-y(j)) + beta*v3(j);
    end
    v2(i) = max(rhs);                           %find maximum
end

%Given v2, solve for v1 
for i=1:length(x)
    for j=1:i
        rhs2(j) = log(x(i)-y(j)) + beta*v2(j);
    end
    v1(i) = max(rhs2);
end










