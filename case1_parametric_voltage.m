% ===========================================================
% PARAMETRIC ESTIMATION: VOLTAGE SIGNAL
% ===========================================================
% This scripts carries out an autoregressive (AR) parametric estimation
% using a PMU-recorded voltage signal for forced oscillation
% identification.
%
% Last modification: 02/07/2023 by SADR.
clc
clear 
close all

% Loading data:
% - `Freq`:     frequency measurements.
% - `fs`:       PMU data sampling rate.
% - `ts`:       PMU data samping time.
% - `time`:     time measurement vector.
% - `t_event`:  start time of the oscillation
% - `Im`:       current magnitude.
% - `Vm`:       voltage magnitude.
load data/sysid_power_case1.mat

%% Step 1: Data Preprocessing

% ---- Extracting the number of samples
N = size(Vm, 1);

% ---- Getting signal (bus 5)
y1 = Vm(1202:1802, 5);
N1 = size(y1, 1);

% ---- Creating the time vector (using the sampling time)
t1 = ts:ts:length(y1)*ts;

% ---- Getting signal (bus 13)
y2 = Vm(1202:1802, 13);
N2 = size(y2, 1);
t2 = ts:ts:length(y2)*ts;

% ---- Detrending measurements
y1 = detrend(y1);
y2 = detrend(y2);

% ---- Taking the number of measurements as N+1 (using N)
N1 = N1 - 1;
N2 = N2 - 1;

% ---- Measurements to from k=1 to k=N
y1_meas = Vm(1203:1802, 5);
y2_meas = Vm(1203:1802, 13);

% ---- Detrending measurements
y1_meas = detrend(y1_meas);
y2_meas = detrend(y2_meas);

%% Step 2: Fitting AR Model

% ---- Minimum and maximum orders of the AR models
na_min = 1;
na_max = 40;

% ---- Initializing containers
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

% -----------------------------------------------------------
% AR estimation
% -----------------------------------------------------------
for na=na_min:na_max
    % ---- Direct solution of the LS estimation problem
    
    % ---- Initialization
    phi1 = zeros(N1, na);
    phi2 = zeros(N2, na);

    % ---- Populating the columns
    for j=1:na
        phi1(j:end-1, j) = -y1(1:end-1-j);
        phi2(j:end-1, j) = -y2(1:end-1-j);
    end
    
    % ---- Regression
    theta1{na} = pinv(phi1)*y1_meas;
    theta2{na} = pinv(phi2)*y2_meas;
    
    % ---- Computation of MSE
    mse1{na} = (1/N1)*(y1_meas - phi1*theta1{na}).'*(y1_meas - phi1*theta1{na});
    mse2{na} = (1/N2)*(y2_meas - phi2*theta1{na}).'*(y2_meas - phi2*theta2{na});

    res1{na} = y1_meas - phi1*theta1{na};
    res2{na} = y2_meas - phi2*theta2{na};

    % ---- MATLAB's method
    sys_AR1{na} = ar(y1, na, 'ls');
    sys_AR2{na} = ar(y2, na, 'ls');
    theta_AR1{na} = sys_AR1{na}.Report.Parameters.ParVector;
    theta_AR2{na} = sys_AR2{na}.Report.Parameters.ParVector;
    res_AR1{na} = y1_meas - phi1*theta_AR1{na};
    res_AR2{na} = y2_meas - phi2*theta_AR2{na};
    mse_AR1{na} = sys_AR1{na}.Report.Fit.MSE;
    mse_AR2{na} = sys_AR2{na}.Report.Fit.MSE;
end

%% Step 3: Model Order Selection

clc
close all

% ---- Plot settings for the MSE error (as a function of the model order)
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

% ---- Plotting the difference between parameters
% ---- of the direct solution and of MATLAB's AR function
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

exportgraphics(fig, 'results/fig_case1_Vm_model-order.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Step 4: Residual Whiteness Analysis

close all

% ---- Settings for the autocorrelation plot
width = 10;
height = 6;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')

% ---- Model order (can be changed according to the results of Step 3))
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

exportgraphics(fig, 'results/fig_case1_Vm_residual-20.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Step 5: Mode frequency computation

close all

% ---- Getting the poles as the roots of the characteristic polynomial
poles_LS = roots([1 theta1{n_a}.']);
poles_LS_CT = log(poles_LS)*fs;
% ---- Converting to CT domain at the sampling frequency
poles_AR = roots([1 theta_AR1{n_a}.']);
poles_AR_CT = log(poles_LS)*fs;

% ---- Settings for the pole plots
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

exportgraphics(fig, 'results/fig_case1_Vm_pole-diagram.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Step 6: Estimation of frequency and damping

% -----------------------------------------------------------
% LSE implementation
% -----------------------------------------------------------

% ---- Computation of the frequency using the resulting poles
w_estimated_LS = imag(poles_LS_CT);

% ---- Computation of the damping using the resulting poles
sigma_estimated_LS = real(poles_LS_CT);
magn_LS = sqrt(w_estimated_LS.^2 + sigma_estimated_LS.^2);

f_estimated_LS = w_estimated_LS / (2*pi)
damp_estimated_LS = (-sigma_estimated_LS ./ magn_LS) * 100

% -----------------------------------------------------------
% AR model
% -----------------------------------------------------------
% ---- Computation of the frequency using the resulting poles
w_estimated_AR = imag(poles_AR_CT);

% ---- Computation of the damping using the resulting poles
sigma_estimated_AR = real(poles_AR_CT);
magn_AR = sqrt(w_estimated_AR.^2 + sigma_estimated_AR.^2);

f_estimated_AR = w_estimated_AR / (2*pi)
damp_estimated_AR = (-sigma_estimated_AR ./ magn_LS) * 100