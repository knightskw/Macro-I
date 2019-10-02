% Shai Knight-Winnig 2016
clear;

% Given from the Problem
Beta=(1/1.04)^0.25;
fi1=0.4515;
fi2=1.4916;
nu1=0.25;
nu2=0.95;
Alpha=0.36;
Delta=0.00788;

% Setting up plots
nk=401; % Capital
nl=101; % Labor
kx=linspace(33,38,nk);
ky=linspace(33,38,nk);
l=linspace(0.3,0.35,nl);

J=100; % Iterations

% Value Function

ll=l; % ll = 1*nl
l=l';
l=l(:,ones(1,nk));
ky=ky(ones(nl,1),:);

% Creating Plot 1
temp_k1=zeros(nk,nk); 
temp_k2=zeros(nk,nk);
temp_k3=zeros(nk,nk);
v1_rhs1=zeros(1,nk);
v1_rhs2=zeros(1,nk);
v1_rhs3=zeros(1,nk);
v1_lhs1=zeros(nk,J);
v1_lhs2=zeros(nk,J);
v1_lhs3=zeros(nk,J);
place_l1=ones(nk,nk);
place_l2=ones(nk,nk);
place_l3=ones(nk,nk);
place_k1=ones(1,nk);
place_k2=ones(1,nk);
place_k3=ones(1,nk);

% Creating Plot 2
ttemp_k1=zeros(nk,nk); 
ttemp_k2=zeros(nk,nk);
ttemp_k3=zeros(nk,nk);
vv1_rhs1=zeros(1,nk);
vv1_rhs2=zeros(1,nk);
vv1_rhs3=zeros(1,nk);
vv1_lhs1=zeros(nk,J);
vv1_lhs2=zeros(nk,J);
vv1_lhs3=zeros(nk,J);
pplace_l1=ones(nk,nk);
pplace_l2=ones(nk,nk);
pplace_l3=ones(nk,nk);
pplace_k1=ones(1,nk);
pplace_k2=ones(1,nk);
pplace_k3=ones(1,nk);

for j=1:J
    v1_rhs1=v1_rhs1(ones(nl,1),:);
    v1_rhs2=v1_rhs2(ones(nl,1),:);
    v1_rhs3=v1_rhs3(ones(nl,1),:);
    vv1_rhs1=vv1_rhs1(ones(nl,1),:);
    vv1_rhs2=vv1_rhs2(ones(nl,1),:);
    vv1_rhs3=vv1_rhs3(ones(nl,1),:);
    
    for i=1:nk
        [temp_k1(i,:),place_l1(i,:)]=max((log(max(((exp(-0.007))*(kx(i)^Alpha).*(l.^(1-Alpha))+...
           +kx(i)*(1-Delta)-ky),1e-12))+...
           fi1*(1/(1-1/nu1))*((1-l).^(1-1/nu1))+Beta*(1/3)*(v1_rhs1+v1_rhs2+v1_rhs3)),[],1);
        [temp_k2(i,:),place_l2(i,:)]=max((log(max((kx(i)^Alpha).*(l.^(1-Alpha))+...
           +kx(i)*(1-Delta)-ky,1e-12))+...
           fi1*(1/(1-1/nu1))*((1-l).^(1-1/nu1))+Beta*(1/3)*(v1_rhs1+v1_rhs2+v1_rhs3)),[],1);
        [temp_k3(i,:),place_l3(i,:)]=max((log(max(((exp(0.007))*(kx(i)^Alpha).*(l.^(1-Alpha))+...
           +kx(i)*(1-Delta)-ky),1e-12))+...
           fi1*(1/(1-1/nu1))*((1-l).^(1-1/nu1))+Beta*(1/3)*(v1_rhs1+v1_rhs2+v1_rhs3)),[],1);
        [v1_lhs1(i,j),place_k1(i)]=max(temp_k1(i,:),[],2);
        [v1_lhs2(i,j),place_k2(i)]=max(temp_k2(i,:),[],2);
        [v1_lhs3(i,j),place_k3(i)]=max(temp_k3(i,:),[],2);
        
        [ttemp_k1(i,:),pplace_l1(i,:)]=max((log(max(((exp(-0.007))*(kx(i)^Alpha).*(l.^(1-Alpha))+...
           +kx(i)*(1-Delta)-ky),1e-12))+...
           fi2*(1/(1-1/nu2))*((1-l).^(1-1/nu2))+Beta*(1/3)*(vv1_rhs1+vv1_rhs2+vv1_rhs3)),[],1);
        [ttemp_k2(i,:),pplace_l2(i,:)]=max((log(max((kx(i)^Alpha).*(l.^(1-Alpha))+...
           +kx(i)*(1-Delta)-ky,1e-12))+...
           fi2*(1/(1-1/nu2))*((1-l).^(1-1/nu2))+Beta*(1/3)*(vv1_rhs1+vv1_rhs2+vv1_rhs3)),[],1);
        [ttemp_k3(i,:),pplace_l3(i,:)]=max((log(max(((exp(0.007))*(kx(i)^Alpha).*(l.^(1-Alpha))+...
           +kx(i)*(1-Delta)-ky),1e-12))+...
           fi2*(1/(1-1/nu2))*((1-l).^(1-1/nu2))+Beta*(1/3)*(vv1_rhs1+vv1_rhs2+vv1_rhs3)),[],1);
        [vv1_lhs1(i,j),pplace_k1(i)]=max(ttemp_k1(i,:),[],2);
        [vv1_lhs2(i,j),pplace_k2(i)]=max(ttemp_k2(i,:),[],2);
        [vv1_lhs3(i,j),pplace_k3(i)]=max(ttemp_k3(i,:),[],2);       
    end
    
    v1_rhs1=v1_lhs1(:,j)';
    v1_rhs2=v1_lhs2(:,j)';
    v1_rhs3=v1_lhs3(:,j)';
    vv1_rhs1=vv1_lhs1(:,j)';
    vv1_rhs2=vv1_lhs2(:,j)';
    vv1_rhs3=vv1_lhs3(:,j)';
end

T=100;

% Random z;
R=rand(1,T);
Z1=(R<=1/3);
Z2=(R>1/3)&(R<=2/3);
Z3=(R>2/3);
Z=1*Z1+2*Z2+3*Z3;

capital1=place_k1;
capital2=place_k2;
capital3=place_k3;
ccapital1=pplace_k1;
ccapital2=pplace_k2;
ccapital3=pplace_k3;

labor1=ones(1,nk);
labor2=ones(1,nk);
labor3=ones(1,nk);
llabor1=ones(1,nk);
llabor2=ones(1,nk);
llabor3=ones(1,nk);

% Main
for i=1:nk
   labor1(i)=place_l1(i,place_k1(i));
   labor2(i)=place_l2(i,place_k2(i));
   labor3(i)=place_l3(i,place_k3(i));
   llabor1(i)=pplace_l1(i,pplace_k1(i));
   llabor2(i)=pplace_l2(i,pplace_k2(i));
   llabor3(i)=pplace_l3(i,pplace_k3(i));
end

Z(1)=2;
optk=ones(1,T);
optl=ones(1,T);
optk(1)=241;
optl(1)=67;

ooptk=ones(1,T);
ooptl=ones(1,T);
ooptk(1)=241;
ooptl(1)=67;


for t=2:T
   if Z(t)==1
       optk(t)=capital1(optk(t-1));
       optl(t)=labor1(optk(t-1));
       ooptk(t)=ccapital1(ooptk(t-1));
       ooptl(t)=llabor1(ooptk(t-1));
   elseif Z(t)==2
       optk(t)=capital2(optk(t-1));
       optl(t)=labor2(optk(t-1));
       ooptk(t)=ccapital2(ooptk(t-1));
       ooptl(t)=llabor2(ooptk(t-1));
   elseif Z(t)==3
       optk(t)=capital3(optk(t-1));
       optl(t)=labor3(optk(t-1));
       ooptk(t)=ccapital3(ooptk(t-1));
       ooptl(t)=llabor3(ooptk(t-1));
   end
end

optk_final=kx(optk);
optl_final=ll(optl);
opty_final=(optk_final.^(Alpha)).*(optl_final.^(1-Alpha));

ooptk_final=kx(ooptk);
ooptl_final=ll(ooptl);
oopty_final=(ooptk_final.^(Alpha)).*(ooptl_final.^(1-Alpha));

% Plotting the Final Graphs
tt=linspace(1,100,100);
figure;
plot(tt,opty_final,'g',tt,oopty_final,'r');
legend('Model1','Model2');
title('y path');

figure;
plot(tt,optl_final,'g',tt,ooptl_final,'r');
legend('Model1','Model2');
title('l path');

