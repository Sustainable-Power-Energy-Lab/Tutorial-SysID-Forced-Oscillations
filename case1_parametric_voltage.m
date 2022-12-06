%% Loading Data
clc
clear 
close all

% Loading given data (measurements in `y`; `N` observations)
load data/sysid_power_case1.mat

N = size(Vm, 1);

%% Loading Data for First Fitting

y1 = Vm(1202:1802, 5);
N1 = size(y1, 1);
t1 = ts:ts:length(y1)*ts;

y2 = Vm(1202:1802,13);
N2 = size(y2, 1);
t2 = ts:ts:length(y2)*ts;

% Detrending measurements
y1 = detrend(y1);
y2 = detrend(y2);

% Taking the number of measurements as N+1 (using N)
% see report for more details
N1 = N1 - 1;
N2 = N2 - 1;

% Measurements to from k=1 to k=N
y1_meas = Vm(1203:1802, 5);
y2_meas = Vm(1203:1802, 13);

y1_meas = detrend(y1_meas);
y2_meas = detrend(y2_meas);

%% Fitting AR Model

% Minimum and maximum orders of the AR models
na_min = 1;
na_max = 40;

% Initializing containers
theta1 = {};
theta2 = {};
mse1 = {};
mse2 = {};
sys_AR1 = {};
theta_AR1 = {};
mse_AR1 = {};
sys_AR2 = {};
theta_AR2 = {};
mse_AR2 = {};

% ================================
% Estimation of the AR model
% ================================
for na=na_min:na_max
    % Building the `phi` matrix
    
    % Initialization
    phi1 = zeros(N1, na);
    phi2 = zeros(N2, na);

    % Populating the columns
    for j=1:na
        phi1(j:end-1, j) = -y1(1:end-1-j);
        phi2(j:end-1, j) = -y2(1:end-1-j);
    end
    
    % Regression
    theta1{na} = pinv(phi1)*y1_meas;
    theta2{na} = pinv(phi2)*y2_meas;
    
    % Computation of MSE
    mse1{na} = (1/N1)*(y1_meas - phi1*theta1{na}).'*(y1_meas - phi1*theta1{na});
    mse2{na} = (1/N2)*(y2_meas - phi2*theta1{na}).'*(y2_meas - phi2*theta2{na});

    res1{na} = y1_meas - phi1*theta1{na};
    res2{na} = y2_meas - phi2*theta2{na};

    % MATLAB's method
    sys_AR1{na} = ar(y1, na, 'ls');
    sys_AR2{na} = ar(y2, na, 'ls');
    theta_AR1{na} = sys_AR1{na}.Report.Parameters.ParVector;
    theta_AR2{na} = sys_AR2{na}.Report.Parameters.ParVector;
    res_AR1{na} = y1_meas - phi1*theta_AR1{na};
    res_AR2{na} = y2_meas - phi2*theta_AR2{na};
    mse_AR1{na} = sys_AR1{na}.Report.Fit.MSE;
    mse_AR2{na} = sys_AR2{na}.Report.Fit.MSE;
end

%% Model Order Selection

clc
close all

width = 10;
height = 6;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')

order_model = na_min:na_max;

subplot(121)
plot(order_model, cell2mat(mse1), 'DisplayName', 'LSE AR', 'color', '#191970')
hold on
plot(order_model, cell2mat(mse_AR1), 'DisplayName', 'MATLAB AR', ...
    'color', 'red', 'linestyle', '--')
title('Mean Squared Error (MSE)');
lgd = legend('interpreter', 'latex');
% lgd.NumColumns = 4;
lgd.Location = 'northeast';
xlim([na_min, na_max])
xlabel('Model Order ($n_a$)', 'interpreter', 'latex');

subplot(122)

diff_par = {};

for i=na_min:na_max
    diff_par{i} = norm(theta1{i} - theta_AR1{i}, 2);
end

plot(na_min:na_max, cell2mat(diff_par), 'color', '#191970')
title('Difference between Parameters')
xlabel('Model Order ($n_a$)', 'interpreter', 'latex');
ylabel('$\Vert\hat{\theta}_{AR}-\hat{\theta}_{ML}\Vert_2$', ...
    'interpreter', 'latex')
xlim([na_min, na_max])

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case1_Vm_model-order.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Residual Whiteness Analysis

close all
width = 10;
height = 6;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')

% Model order
n_a = 20;

subplot(221)
plot(linspace(1/fs,N1/fs,N1), res1{n_a}, 'color', '#191970')
title_string = sprintf('LSE Residual ($n_a$ = %d)',...
    n_a);

title(title_string, 'interpreter', 'latex')
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$r[t]$', 'interpreter', 'latex')

subplot(223)
autocorr(res1{n_a}, round(N1/4))

subplot(222)
plot(linspace(1/fs,N1/fs,N1), res_AR1{n_a}, 'color', '#191970')
title_string = sprintf('AR Residual ($n_a$ = %d)',...
    n_a);
title(title_string, 'interpreter', 'latex')
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$r[t]$', 'interpreter', 'latex')

subplot(224)
autocorr(res_AR1{n_a}, round(N1/4))


exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case1_Vm_residual-20.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Mode frequency computation

close all

poles_LS = roots([1 theta1{n_a}.']);
poles_LS_CT = log(poles_LS)*fs;

poles_AR = roots([1 theta_AR1{n_a}.']);
poles_AR_CT = log(poles_LS)*fs;

width = 8;
height = 6;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')
scatter(real(poles_LS_CT), imag(poles_LS_CT), 100, 'x', ...
    'DisplayName','LS-based', 'color', '#191970')
hold on
scatter(real(poles_AR_CT), imag(poles_AR_CT), 100, 'x', 'color', 'red', ...
    'DisplayName', 'AR')

lgd = legend('interpreter', 'latex');
% lgd.NumColumns = 4;
lgd.Location = 'northwest';
title('Pole diagram')
xlabel('Real Axis ($\sigma$)', 'interpreter', 'latex')
ylabel('Imaginary Axis ($j\omega$)', 'interpreter', 'latex')
xlim([-20 0.1])

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case1_Vm_pole-diagram.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Estimation of frequency and damping

w_estimated_LS = imag(poles_LS_CT);
sigma_estimated_LS = real(poles_LS_CT);

magn_LS = sqrt(w_estimated_LS.^2 + sigma_estimated_LS.^2);

f_estimated_LS = w_estimated_LS / (2*pi)
damp_estimated_LS = (-sigma_estimated_LS ./ magn_LS) * 100

w_estimated_AR = imag(poles_AR_CT);
sigma_estimated_AR = real(poles_AR_CT);

magn_AR = sqrt(w_estimated_AR.^2 + sigma_estimated_AR.^2);

f_estimated_AR = w_estimated_AR / (2*pi)
damp_estimated_AR = (-sigma_estimated_AR ./ magn_LS) * 100
