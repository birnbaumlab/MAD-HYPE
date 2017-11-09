function holec_patrick_huang

%% NOTE:
% I learned later on that k5,k6 are defined as 1.7*k1, so all simulations
% still technically solver for this variable, but it is disconnected from
% the ODE (therefore, reported values are meaningless

%% Variable setup

close all;

% data estimates
K1 = 0.0981; k2 = 0.0389; k3 = 0.0387; k4 = 0.0043;

%kb12 = 0.0110; kb2 = 0.025; kb3 = 0.0105; % original guess
kb12 = .0099; kb2 = .0075; kb3 = 1.89e-8; % modified guess (makes faster)

% unknowns
K5 = 0.1822; k6 = 0.0758;

% data pulling
data = csvread('plasma_totalF18_FIG3.csv',1,0);
t_ra = data(:,1); y_ra = data(:,2);
data = csvread('plasma_3OMFDFraction_FIG4.csv',1,0);
t_data = data(:,1); y_data = data(:,2);


%% Fig. 4

figure('units','normalized','outerposition',[0 0 .6 .5]);

% fits parameters kb12,kb2,kb3
disp('Starting parameter fitting for Fig. 4...');
params = nlinfit(t_data,y_data,@nlinfit1,[kb12,kb2,kb3]);

% redo ODE with best fit parameters
options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
p = {params(1),params(2),params(3),t_ra,y_ra};
[~,y] = ode15s(@odefun1,t_data,[0,0],options,p);

disp(['kb12: ',num2str(p{1})]);
disp(['kb2: ',num2str(p{2})]);
disp(['kb3: ',num2str(p{3})]);

% find normalizing vector for F18 signal
y_norm = interp1(t_ra,y_ra,t_data);

% normalize fraction of 3-OMFD
y_omfd_frac = y(:,1)./y_norm; y_omfd_frac(1) = 0;

% plot
plot(t_data,y_omfd_frac); hold on;
plot(t_data,y_data,'sr');

% plot formatting
xlabel('Time (min)');
ylabel('3-OMFD fraction');
legend({'Fitting result','Measured data'},'Location','southeast');

title(['\bf{Figure 4: }' ...
    '\rm{3-OMFD fraction of F18 signal after FDOPA injection. Experimental results were fit to the' char(10) ...
    'model described by Huang et al, producing constants for k_{b12},k_{b2}, and k_{b3}. Fit is almost' char(10) ...
    'exactly in line with input dataset.' char(10) ...
    '}']);


%% Fig. 5

figure('units','normalized','outerposition',[0 0 .6 .5]);

% pull some data
data = csvread('plasma_totalF18_FIG5.csv',1,0);
t_5 = data(:,1); y_5 = data(:,2);

% set up passable variable of constants
p = {kb12,kb2,kb3,t_5,y_5};

% complete ODE to solve for concentraiton of species
options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
[~,y] = ode15s(@odefun1,t_5,[0,0],options,p);

% calculate fraction of each
y_omfd_frac = y(:,1)./y_5;

% switch to concentraiton of species (yes, I know this is redundant)
FDOPA = (1 - y_omfd_frac).*y_5;
OMFD = y_omfd_frac.*y_5;

% plot and format
semilogy(t_5,y_5,'s-r'); hold on;
semilogy(t_5,FDOPA,'--k');
semilogy(t_5,OMFD,':b');

xlabel('Time (min)');
ylabel('Radioactivity (uCi/ml)');
legend({'Total F-18','FDOPA','3-OMFD'});
xlim([0 120]); ylim([0.001 10]);

title(['\bf{Figure 5: }' ...
    '\rm{Radioactivity resulting from FDOPA and 3-OMFD after F18 injection. Results from Fig. 4 fraction' char(10) ...
    'were extrapolated to new set of plasma F18 signal. Intially, F18 signal is dominated by FDOPA but' char(10) ...
    'switches dominant fraction around 30 minutes.' char(10) ...
    '}']);

%% Fig. 6

figure('units','normalized','outerposition',[0 0 .6 .5]);

% Extracting data from files
data = csvread('tissue_radioactivity_cerebellum_FIG6.csv',1,0);
t_f18_cb = data(:,1); y_f18_cb = data(:,2);

data = csvread('tissue_radioactivity_striatum_FIG6.csv',1,0);
t_f18_st = data(:,1); y_f18_st = data(:,2);

data = csvread('plasma_totalF18_FIG5.csv',1,0);
t_5_plasma = data(:,1); y_5_plasma = data(:,2);

data = csvread('plasma_totalF18_FIG3.csv',1,0);
t_3_plasma = data(:,1); y_3_plasma = data(:,2);

data = csvread('wholeblood_totalF18_FIG3.csv',1,0);
t_3_wb = data(:,1); y_3_wb = data(:,2);


% Setting the table for p1/p2 fitting
p = {kb12,kb2,kb3,t_3_plasma,y_3_plasma};

% solve OD for fraction species
options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
[~,y] = ode15s(@odefun1,t_3_plasma,[0,0],options,p);

y_omfd_frac = y(:,1)./y_3_plasma;

% resulting concentrations of FDOPA and 3-OMFD
FDOPA = (1 - y_omfd_frac).*y_3_plasma;
OMFD = y_omfd_frac.*y_3_plasma;

% set up parameter package
p_extra = {t_3_plasma,FDOPA,OMFD};

% fit for p1,p2 parameters in whole blood
disp('Starting parameter fitting for p1,p2 in Fig. 3...)');
params = nlinfit(t_3_wb,y_3_wb,@(p,t)nlinfit2(p,t,p_extra),[0.7,0.5]);
p1 = params(1); p2 = params(2);
disp('Finished fitting.');

disp(['p1: ',num2str(p1)]);
disp(['p2: ',num2str(p2)]);

% Setting the table for FDOPA/3-OMFD concentration
p = {kb12,kb2,kb3,t_5_plasma,y_5_plasma};

% resolve system with p1,p2 values in mind, on fig. 5 data now
options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
[~,y] = ode15s(@odefun1,t_5_plasma,[0,0],options,p);

% calculate fraction
y_omfd_frac = y(:,1)./y_5_plasma;

% calculate concentration of fig. 5 data
FDOPA = (1 - y_omfd_frac).*y_5_plasma/p1;
OMFD = y_omfd_frac.*y_5_plasma/p2;

% Final fitting

% Striatum analysis
p_extra = {p1,p2,t_5_plasma,FDOPA,OMFD};

% fit K1,k2,k3,K5,k6 parameters
disp('Starting parameter fitting for striatum model in Fig. 6B...');
params = nlinfit(t_f18_st,y_f18_st,@(p,t)nlinfit3(p,t,p_extra),[K1,k2,k3,K5,k6]);
disp('Finished fitting.');

% solve differential with known constants
p0 = {params(1),params(2),params(3),params(4),params(5),t_5_plasma,FDOPA,OMFD};
[t,y] = ode15s(@odefun2,t_f18_st,[0,0,0],options,p0);

% calculate total concentration of species
C_wb = interp1(t,y(:,1),t) + interp1(t,y(:,3),t) + interp1(t,y(:,2),t);

% plot radioactivity of particular species
plot(t_f18_st,y_f18_st,'sr'); hold on;
plot(t_f18_st,C_wb,'-r');
plot(t_f18_st,p2*y(:,3),'--k');
plot(t_f18_st,p1*y(:,1) + y(:,3),':k');

disp(['K1: ',num2str(params(1))]);
disp(['K2: ',num2str(params(2))]);
disp(['K3: ',num2str(params(3))]);
disp(['K5: ',num2str(params(4))]);
disp(['K6: ',num2str(params(5))]);

% Cerebellum analysis
p_extra = {p1,p2,t_5_plasma,FDOPA,OMFD};

% fit K1,k2,k3,K5,k6 parameters
disp('Starting parameter fitting for cerebellum model in Fig. 6B...');
params = nlinfit(t_f18_cb,y_f18_cb,@(p,t)nlinfit4(p,t,p_extra),[K1,k2,K5,k6]);
disp('Finished fitting.');

% solve differential with known constants
p0 = {params(1),params(2),params(3),params(4),t_5_plasma,FDOPA,OMFD};
[t,y] = ode15s(@odefun3,t_f18_cb,[0,0],options,p0);

% calculate total concentration of species
C_wb = p1*0.6*interp1(t,y(:,1),t) + p2*interp1(t,y(:,2),t);

% plot radioactivity of particular species
plot(t_f18_cb,y_f18_cb,'ob'); hold on;
plot(t_f18_cb,C_wb,'-b');

disp(['Kp1: ',num2str(params(1))]);
disp(['Kp2: ',num2str(params(2))]);
disp(['Kp5: ',num2str(params(3))]);
disp(['Kp6: ',num2str(params(4))]);

% Plot Labelling
xlabel('Time (min)');
ylabel('Radioactivity (uCi/ml)');
legend({'Striatum (data)','Striatum (model)','Striatum (3-OMFD)','Striatum (3-OMFD + FDOPA)','Cerebellum (data)','Cerebellum (model)'});

title(['\bf{Figure 6B: }' ...
    '\rm{Concentration of various species after time course injection in cerebellum and striatum. The' char(10) ...
    'concentration of FDA in the striatum can be interpreted as the difference between the 3-OMFD + FDOPA' char(10) ...
    'signal and the total signal in the striatum.' char(10) ...
    '}']);

%% Fig. 7 (Modification)

% We are going to be looking at the systems resiliency to noise being added
% to the system input

figure('units','normalized','outerposition',[0 0 .6 1]);
noise = [0,0.01,0.02,0.05,0.1];
trials = 25;

% This code takes literally hours to run, so unless you have patience,
% trust my respresentative sample is legit below
if 0
    k1_all = zeros(size(noise,2),trials);
    k2_all = zeros(size(noise,2),trials);
    k3_all = zeros(size(noise,2),trials);

    % fits parameters kb12,kb2,kb3
    disp('Starting parameter fitting series for Fig. 7...');
    for n = 1:length(noise)
        disp(['Starting trial ',num2str(n)])
        for t = 1:trials
            params = nlinfit(t_data,abs(y_data + normrnd(0,y_data*noise(n),size(y_data))),@nlinfit1,[kb12,kb2,kb3]);
            k1_all(n,t) = params(1); k2_all(n,t) = params(2); k3_all(n,t) = params(3); 
        end
    end


else % else if you dont want to wait forever...
    k1_all = [0.00990188255976311,0.00990188255976311,0.00990188255976311;0.0101354506498536,0.00990573134932737,0.0100749123149245;0.0100705232187096,0.00961585563411375,0.0101658141836586;0.00984074630029715,0.00967705287781844,0.00990088913680848;0.0122197626749377,0.00961581463689033,0.0100387677608780];
    k2_all = [0.00750016933242117,0.00750016933242117,0.00750016933242117;0.00795825053428397,0.00746349355351963,0.00800058129222460;0.00794404706358343,0.00723144404458605,0.00773431927929209;0.00801900174395646,0.00738043349298682,0.00806745873642345;0.00795983948301586,0.00974801383663627,0.0109903779975501];
    k3_all = [6.90821039520759e-08,6.90821039520759e-08,6.90821039520759e-08;4.80599101640014e-08,-1.11743143415402e-07,-6.60765090643153e-08;-3.17509691636442e-08,7.40687983815870e-08,-5.30843673554981e-08;-5.58815884311175e-08,-0.000540484669425541,5.71533312345032e-08;-1.06241375837545e-07,0.00654082261103412,0.00709303878257514];
end

k = {abs(k1_all),abs(k2_all),abs(k3_all)};

for i = 1:length(k)
    subplot(3,1,i);
    for j = 2:length(noise)
        bar(mean(k{i},2),'FaceColor',[0 .5 .5],'EdgeColor',[0 0 0],'LineWidth',2); hold on;
    end
    %set(gca,'yscale','log')
    set(gca,'XTickLabel',{'0% noise', '1% noise', '2% noise', '5% noise', '10% noise'});
    disp(std(k{i},0,2))
    errorbar(1:5,mean(k{i},2),std(k{i},0,2),'.k')
end

subplot(3,1,1);
ylim([0.005 0.015]); ylabel('K_{1} value (ml/min/g)');
subplot(3,1,2);
ylabel('k_{2} value (/min)');
subplot(3,1,3);
ylim([0 6e-3]); ylabel('k_{3} value (/min)');
subplot(3,1,1);
title(['\bf{Figure 7: }' ...
    '\rm{Sensitivity of FDOPA/3-OMFD conversion to noise. As noise is added to the experimental data,' char(10) ...
    'k_{1} and k_{2} remain largely resistant to perturbation in average rate value. However, k_{3}' char(10) ...
    'radically shifts with slightest addition of noise suggesting this may be an extraneous degree' char(10) ...
    'of freedom.' char(10) ...
    '}']);

%% Fig. 8

figure('units','normalized','outerposition',[0 0 .6 .5]);

% fits parameters kb12,kb2,kb3
disp('Starting parameter fitting for Fig. 4...');
params = nlinfit(t_data,y_data,@nlinfit5,[kb12,kb2,kb3]);

% redo ODE with best fit parameters
options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
p = {params(1),params(2),params(3),t_ra,y_ra};
[~,y] = ode15s(@odefun1,t_data,[0,0],options,p);

disp(['kb12: ',num2str(p{1})]);
disp(['kb2: ',num2str(p{2})]);
disp(['kb3: ',num2str(p{3})]);

% find normalizing vector for F18 signal
y_norm = interp1(t_ra,y_ra,t_data);

% normalize fraction of 3-OMFD
y_omfd_frac = y(:,1)./y_norm; y_omfd_frac(1) = 0;

% plot
plot(t_data,y_omfd_frac); hold on;
plot(t_data,y_data,'sr');

% plot formatting
xlabel('Time (min)');
ylabel('3-OMFD fraction');
legend({'Fitting result','Measured data'},'Location','southeast');

title(['\bf{Figure 8: }' ...
    '\rm{Fitting of 3-OMFD fraction when k_{b12} constant is omitted.' char(10) ...
    '' char(10) ...
    '}']);
end




%% Nonlinear fitting solvers

% This is for solving for the kb12, kb2, kb3 values for Fig. 4/5
function y = nlinfit1(p,t)
% dont read a file since its faster to declare here
t_ra = [ 0;    0.1095;    0.2987;    0.5256;    0.7586;    0.9834;    1.2080;    1.4174;    1.6831;    1.9300;    2.1637;    2.9453;    4.1013;    5.0370;    7.2585;   10.1033;   20.2575;   30.1943;   61.9289;   94.4310;    121.3510];
y_ra = [ 0;    0.0001;    0.0068;    1.9278;    1.6733;    0.6337;    0.5741;    0.5343;    0.4396;    0.4248;    0.3900;    0.3155;    0.2561;    0.2164;    0.1470;    0.1325;    0.0937;    0.0743;    0.0505;    0.0362;    0.0316];

y_data = interp1(t_ra,y_ra,t);
options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
p = {abs(p(1)),abs(p(2)),abs(p(3)),t_ra,y_ra};
[~,y] = ode15s(@odefun1,t,[0,0],options,p);
y = y(:,1)./y_data;
end

% This is for solving for the kb12, kb2, kb3 values for Fig. 4/5
function C_wb = nlinfit2(p,t,p_extra)

t_data = p_extra{1};
FDOPA = p_extra{2};
OMFD = p_extra{3};
C_wb = p(1)*0.6*interp1(t_data,FDOPA,t) + p(2)*interp1(t_data,OMFD,t);

end

% This is for solving for the striatum constants
function C_wb = nlinfit3(p,t,p_extra)
p1 = p_extra{1}; p2 = p_extra{2};
t_data = p_extra{3}; FDOPA = p_extra{4}; OMFD = p_extra{5};
p0 = {abs(p(1)),abs(p(2)),abs(p(3)),abs(p(4)),abs(p(5)),t_data,FDOPA,OMFD}; % pass to function

options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
[~,y] = ode15s(@odefun2,t,[0,0,0],options,p0);
C_wb = interp1(t,y(:,1),t) + interp1(t,y(:,3),t) + interp1(t,y(:,2),t);
end

% This is for solving for the cerebellum constants
function C_wb = nlinfit4(p,t,p_extra)
p1 = p_extra{1}; p2 = p_extra{2};
t_data = p_extra{3}; FDOPA = p_extra{4}; OMFD = p_extra{5};
p0 = {p(1),p(2),p(3),p(4),t_data,FDOPA,OMFD}; % pass to function

options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
[~,y] = ode15s(@odefun3,t,[0,0],options,p0);
C_wb = p1*0.6*interp1(t,y(:,1),t) + p2*interp1(t,y(:,2),t);
end

% This is for solving for the kb12, kb2, kb3 values for Fig. 4/5
function y = nlinfit5(p,t)
% dont read a file since its faster to declare here
t_ra = [ 0;    0.1095;    0.2987;    0.5256;    0.7586;    0.9834;    1.2080;    1.4174;    1.6831;    1.9300;    2.1637;    2.9453;    4.1013;    5.0370;    7.2585;   10.1033;   20.2575;   30.1943;   61.9289;   94.4310;    121.3510];
y_ra = [ 0;    0.0001;    0.0068;    1.9278;    1.6733;    0.6337;    0.5741;    0.5343;    0.4396;    0.4248;    0.3900;    0.3155;    0.2561;    0.2164;    0.1470;    0.1325;    0.0937;    0.0743;    0.0505;    0.0362;    0.0316];

y_data = interp1(t_ra,y_ra,t);
options = odeset('RelTol',1e-13,'AbsTol',[1e-13]);
p = {abs(p(1)),abs(p(2)),0,t_ra,y_ra}; % look! p3 doesnt exist anymore
[~,y] = ode15s(@odefun1,t,[0,0],options,p);
y = y(:,1)./y_data;
end

%% ODE Function Solvers

% This is for solving the C_omfd fraction
function ydot = odefun1(t,y,p)

%paramsCell = num2cell(p);
[kb12,kb2,kb3] = p{1:end-2};
t_ra = p{end-1}; y_ra = p{end};

yCell = num2cell(y);
[C_omfd,C_x] = yCell{:};

C_fd = interp1(t_ra,y_ra,t);

dC_omfd = kb12*C_fd - (kb12 + kb2)*C_omfd + kb3*C_x;
dC_x = kb2*C_omfd - kb3*C_x;

ydot = [dC_omfd;dC_x];

end

% This is for solving the striatum concentraitons (Fig. 6B)
function ydot = odefun2(t,y,p)

[K1,k2,k3,K5,k6] = p{1:5}; 
k4 = 0; % by fig.6b description
t_data = p{6}; FDOPA = p{7}; OMFD = p{8};

yCell = num2cell(y);
[C_1,C_2,C_3] = yCell{:};

C_fd = interp1(t_data,FDOPA,t);
C_omfd = interp1(t_data,OMFD,t);

dC_1 = K1*C_fd - k2*C_1 - k3*C_1;
dC_2 = k3*C_1 - k4*C_2;
dC_3 = 1.7*K1*C_omfd - 1.7*K1*C_3;

ydot = [dC_1,dC_2,dC_3]';

end

% This is for solving the cerebellum concentraitons (Fig. 6B)
function ydot = odefun3(t,y,p)

[K1,k2,K5,k6] = p{1:4}; 
t_data = p{5}; FDOPA = p{6}; OMFD = p{7};

yCell = num2cell(y);
[C_1,C_2] = yCell{:};

C_fd = interp1(t_data,FDOPA,t);
C_omfd = interp1(t_data,OMFD,t);

dC_1 = K1*C_fd - k2*C_1;
dC_2 = 1.7*K1*C_omfd - 1.7*K1*C_2;

ydot = [dC_1,dC_2]';

end



