% Shai Knight-Winnig 2016
% This file computes Arrow securities for two agents with heterogeneus beliefs

%% Parameters

% y individual shocks matrix 2x2 rows are periods and columns states.


% Shocks 
% individual shocks columns are the shocks for each individual 

%%%%%%%%%%%%%%%%%%%%%%%% part I change here for diferent shocks pelase be sure
%%%%%%%%%%%%%%%%%%%%%%%% that all shocks add up to 1 in the columns, change
%%%%%%%%%%%%%%%%%%%%%%%% also beliefs, also the number of grids and limits
%%%%%%%%%%%%%%%%%%%%%%%% of the consumption , notice that the income shock
%%%%%%%%%%%%%%%%%%%%%%%% and the consumption are related, it mean you can
%%%%%%%%%%%%%%%%%%%%%%%% not set up a income shock of 1 and consmp of 10.

%This are the endowments and state up and down in the columns and periods
%in the rows.
y = [
    0.7 0.3
    0.7 0.3
    ];

% Agregatte shocks (This sum all the rows so it ends up as a vector of 1)
Y = sum(y,2);

% Beliefs 
% Matrix of probabilities for
    %P= pessimistic/real agent This believe should
    %O= Optimistic for good states
PiP = [
    0.7 0.3
    0.7 0.3
    ];

% Matrix of probabilities of optimistic agent this should be higher for
% good
PiO = [
    0.8 0.2
    0.8 0.2
    ];

SPts = 2;
Bet = 0.9;
Gam = 1.5;
Gam_i= 1/(Gam)

%%%%%%%%% change the number of points of the consumption
% state space from 0 to 1 21 points
Min_cg= 1e-10;
Max_cg = 1.000;
points_cg = 50; % increase number of points for precision
CGrid = linspace(Min_cg, Max_cg, points_cg);




%% Question a. 
    % Numerically solve for the wealth distribution at(st) and the pricing kernel Qt (s'|st) as functions of: 
    %   -*)  Consumption
    %   -**) Exogenous state st. 
    % Show that there is an one-to-one mapping from the consumption  allocation to the wealth distribution given the exogenous state.




% Q and Cp=C_t+1 have explicit solutions in terms of previous consumption see results from 3e and 3f
% Lets find Q and C_t+1 for each C_t, s_t+1 and S_t 
for i=1:points_cg
    
    %Defining consumption numerically independent of any state
    ci = CGrid(i);
    
    for s=1:SPts
        for sp=1:SPts
            cj=Y(s)-ci
            Q(s,sp,i) = Bet * ( ( (PiP(s,sp)^(Gam_i))*ci/Y(sp)) + ((PiO(s,sp)^(Gam_i))*(cj)/Y(sp)))^Gam;
            Cp(s,sp,i) = Y(sp) * (  (PiP(s,sp)^(Gam_i)*ci) / ((PiP(s,sp)^(Gam_i)*ci + PiO(s,sp)^(Gam_i)*cj))  );
        end
    end
end


% knowing how much arrow security the pessimist have will be sufficient to
% represent wealth distribution. Remember the pesimist is the consumption
% of c_i defined by our grid and in period 0 tehy do not have any arro
% securities

% Asset matrix. 2x21 Matrix (2 states 21 consumption grid) 
A_t_pes = zeros(2,points_cg);

% iterate until we find a convergence between the asset holdings of each period
Metric = 1;
Tolerance = 1e-6;
while (Metric > Tolerance)
    
% for each point in the consumption grid , s and s' we compute the mapping
% from c_t to a given the shock s_t
    for i=1:points_cg
        ci = CGrid(i);
        
         for s=1:SPts
            
            % for each state of period t (first column of s) we compute the
            % a(s_t) is positive if c>y(s_t)
            
            temp_A_t_pes(s,i) = ci-y(s,1);
            for sp=1:SPts
                %prices (2x2x21) along the diferent 21 grids.
                q = Q(s,sp,i);
                %consumption of the pesimist under diferent states and grids (2x2x21)
                cp = Cp(s,sp,i);
                       
                temp_A_t_pes(s,i) = temp_A_t_pes(s,i) + q*interp1(CGrid,A_t_pes(sp,:),cp,'spline','extrap');
            end
        end 
    end
    % Computing the metric for convergence
    Metric = max(abs(temp_A_t_pes(:)-A_t_pes(:)));
    display(Metric);
    
    %Asset holding for each state of period t for diferent grid of consumption (2x21)
    A_t_pes = temp_A_t_pes;

end

% compute current consumption to future asset
for s=1:SPts
    for sp=1:SPts
        % Future asset holding of pesimist 
        A_t1_pess(s,sp,:) = interp1(CGrid,A_t_pes(s,:),Cp(s,sp,:),'spline','extrap');
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%% III part In general you can change the titles legen colors of the graphs
%%%%%%%%%%%%%%%%%%%%%%%% dont touch the squeeze or the hold on and figure
%%%%%%%%%%%%%%%%%%%%%%%% unless have a certificate in matlab , I did it and
%%%%%%%%%%%%%%%%%%%%%%%% screw everything twice :(


%% Question b.

% b) next period arrow security of the pessimist (for s'=U,D) as functions of today's arrow security of the pessimist. 
 
figure(1);
hold on;
plot(squeeze(A_t_pes(1,:)), squeeze(A_t1_pess(1,1,:)),'g');
plot(squeeze(A_t_pes(1,:)), squeeze(A_t1_pess(1,2,:)),'r');
plot(squeeze(A_t_pes(1,:)), squeeze(A_t_pes(1,:)),'k-');
legend('Up','Down','a(S_{t+1})=a(S_{t})');
xlabel('a(S_{t}) - P. agent');
ylabel('a(S_{t+1}) - P. agent');

%% Question c. 

% c)  pricing kernels as functions of the current wealth distribution and the exogenous state.

figure(2);
hold on;
plot(squeeze(A_t_pes(1,:)), squeeze(Q(1,1,:)),'g');
plot(squeeze(A_t_pes(1,:)), squeeze(Q(1,2,:)),'r');
legend('Up','Down');
xlabel('Asset holdings - P agent');
ylabel('Asset Price - P. agent');

%{
%}

%% Question d. 

% generating shocks IID
N = 1000;
rn = rand(1,N);
[~,St] = histc(rn,cumsum(PiP(1,:)));
St = St+1;

% starting point 
Ct(1) = 0.1;

% simulate forward
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