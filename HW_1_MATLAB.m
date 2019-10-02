% Shai Knight-Winnig 2016
% Initialize

%Endowments
y = [
    0.7 0.3
    0.7 0.3
    ];

% Shocks 
Y = sum(y,2);

%P = Pessimist
%O = Optimist
    
Prob_P = [
    0.7 0.3
    0.7 0.3
    ];

% Optimist probabilities

Prob_O = [
    0.8 0.2
    0.8 0.2
    ];

SPts = 2;
Beta = 0.9;
Gamma = 1.5;
Gamma_i= 1/(Gamma)

% State Space
Min_cg= 1e-10;
Max_cg = 1.000;
points_cg = 50; % increase number of points 
CGrid = linspace(Min_cg, Max_cg, points_cg);

% Question a. 

% 3e and 3f

for i=1:points_cg
    
    %Consumption
    cons_i = CGrid(i);
    
    for s=1:SPts
        for sp=1:SPts
            cj=Y(s)-cons_i
            Q(s,sp,i) = Beta * ( ( (Prob_P(s,sp)^(Gamma_i))*cons_i/Y(sp)) + ((Prob_O(s,sp)^(Gamma_i))*(cj)/Y(sp)))^Gamma;
            Cp(s,sp,i) = Y(sp) * (  (Prob_P(s,sp)^(Gamma_i)*cons_i) / ((Prob_P(s,sp)^(Gamma_i)*cons_i + Prob_O(s,sp)^(Gamma_i)*cj))  );
        end
    end
end

% Asset matrix
A_t_pes = zeros(2,points_cg);

% Iterate
Metric = 1;
Tolerance = 1e-6;
while (Metric > Tolerance)
    
% Map schocks s and s'
    for i=1:points_cg
        cons_i = CGrid(i);
        
         for s=1:SPts
            
            % a(s_t)
            temp_A_t_pes(s,i) = cons_i-y(s,1);
            for sp=1:SPts
                % Prices 
                q = Q(s,sp,i);
                % Consumption of the Pesimist
                cp = Cp(s,sp,i);                      
                temp_A_t_pes(s,i) = temp_A_t_pes(s,i) + q*interp1(CGrid,A_t_pes(sp,:),cp,'spline','extrap');
            end
        end 
    end
    % Compute Convergence
    Metric = max(abs(temp_A_t_pes(:)-A_t_pes(:)));
    display(Metric);
    
    %Asset Holding
    A_t_pes = temp_A_t_pes;

end

% Compute Consumption 
for s=1:SPts
    for sp=1:SPts
        % Future Asset
        A_t1_pess(s,sp,:) = interp1(CGrid,A_t_pes(s,:),Cp(s,sp,:),'spline','extrap');
        
    end
end

% Question b.

% b) Next Period Arrow Security
figure(1);
hold on;
plot(squeeze(A_t_pes(1,:)), squeeze(A_t1_pess(1,1,:)),'g');
plot(squeeze(A_t_pes(1,:)), squeeze(A_t1_pess(1,2,:)),'r');
plot(squeeze(A_t_pes(1,:)), squeeze(A_t_pes(1,:)),'k-');
legend('Up','Down','a(S_{t+1})=a(S_{t})');
xlabel('a(S_{t}) - P. agent');
ylabel('a(S_{t+1}) - P. agent');

% Question c 

% c)  Pricing Kernels 
figure(2);
hold on;
plot(squeeze(A_t_pes(1,:)), squeeze(Q(1,1,:)),'g');
plot(squeeze(A_t_pes(1,:)), squeeze(Q(1,2,:)),'r');
legend('Up','Down');
xlabel('Asset holdings - P agent');
ylabel('Asset Price - P. agent');

%{
%}

% Question d. 

% Shocks IID
N = 1000;
rn = rand(1,N);
[~,St] = histc(rn,cumsum(Prob_P(1,:)));
St = St+1;

% Initialize 
Ct(1) = 0.1;

% Iterate
for t=1:N
    if (t<N)
        Ct(t+1) = interp1(CGrid,squeeze(Cp(St(t),St(t+1),:)),Ct(t),'spline','extrap');
    end
    for sp=1:SPts
        Qt(sp,t) = interp1(CGrid,squeeze(Q(St(t),sp,:)),Ct(t),'spline','extrap');
    end
    At(t) = interp1(CGrid,squeeze(A_t_pes(St(t),:)),Ct(t),'spline','extrap');
end

%  Plot wealth distribution of pessimist and time
figure(3);
plot(At,'b');
xlabel('time');
ylabel('a_{t} - P. agent');

% Plot asset price for each state and time.
figure(4);
plot(Qt(1,:),'g');
hold on;
plot(Qt(2,:),'r');
legend('Up','Down');
xlabel('time');
ylabel('Q(s_{t+1}|s_{t})');


%end