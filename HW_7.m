% Shai Knight-Winnig 2016
clear;

%Values Calculated in Question 1
beta=75/77;
A=1;
alpha=3/10;
delta=1/15;
tau_star=0.2;
K_star=4.804;
c_star=1.249;

T=100;  %given

u_prime=1/c_star;
u_double_prime=-c_star^(-2);
F=@(K) K^alpha
W=@(K) (1-alpha)*K^alpha;
F1=alpha*K_star^(alpha-1);
F11=alpha*(alpha-1)*K_star^(alpha-2);
w_prime=alpha*(1-alpha)*K_star^(alpha-1);

syms c1 c2

H1=w_prime+(1+(F1-delta)*(1-tau_star))-c1+K_star*F11*(1-tau_star);
H2=-K_star*(F1-delta)-c2;

eq1=u_double_prime*c1-beta*u_double_prime*(1+(F1-delta)*(1-tau_star))*c1*H1-beta*u_prime*F11*H1*(1-tau_star);
eq2=u_double_prime*c2-beta*u_double_prime*(1+(F1-delta)*(1-tau_star))*(c1*H2+c2)-beta*u_prime*(F11*H2*(1-tau_star)-(F1-delta))

[c1_s, c2_s] = solve([eq1, eq2]);

c1_ss=double(c1_s(2));
c2_ss=double(c2_s(2));

c_bar=@(K,tau) c_star+(c1_ss*(K-K_star)+c2_ss*(tau-tau_star));
H_bar=@(K,tau) W(K)+K*(1+(alpha*K^(alpha-1)-delta)*(1-tau))-c_bar(K,tau);

KK(1)=K_star;
YY(1)=1.601;
WW(1)=W(K_star);
%Plotting {K,Y,W}
for i=2:T
    KK(i)=H_bar(KK(i-1),0.1);
    YY(i)=F(KK(i));
    WW(i)=W(KK(i));
end

figure(2), plot(KK,'Linewidth',3);
figure(2), title('Change in Tax Rate From 0.2 to 0.1');
figure(2), xlabel('t');
figure(2), ylabel('K_t');

xgrid=linspace(0,7);
c_barr=@(K,tau) c_star+(c1_ss.*(K-K_star)+c2_ss.*(tau-tau_star));
H_barr=@(K,tau) (1-alpha).*K.^alpha+K.*(1+(alpha*K.^(alpha-1)-delta).*(1-tau))-c_barr(K,tau);

figure(1), plot(xgrid,H_barr(xgrid,0.2),'Linewidth',1);
figure(1), hold on;
figure(1), plot(xgrid,H_barr(xgrid,0.1),'Linewidth',1);
figure(1), plot(xgrid,xgrid,'Linewidth',1);
figure(1), legend('tax rate=0.2','tax rate=0.1','45 degree line');
figure(1), xlabel('K');
figure(1), ylabel('K Prime');
figure(1), hold off;









