% ===========================================================
% NON-PARAMETRIC ESTIMATION: VOLTAGE SIGNAL
% ===========================================================
% This scripts carries out a Welch-based parametric estimation
% using a PMU-recorded voltage signal for forced oscillation
% identification.
%
% Last modification: 02/06/2023 by SADR.

clc
clear
close all

% ---- Loading given data
load('data\sysid_power_case1.mat');

%% Step 0: Visualization

% ---- Setting up plotting parameters
width = 8;
height = 6;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')

subplot(211)
plot(time, Vm(:, 5), 'color', '#191970')
xlim([0, max(time)])
title('Voltage at Bus 5')
xlabel('Time (s)')
ylabel('$\vert V \vert$ (kV)', 'interpreter', 'latex')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

subplot(212)
plot(time, Vm(:, 13),  'color', '#191970')
xlim([0, max(time)])
title('Voltage at Bus 13')
xlabel('Time (s)')
ylabel('$\vert V \vert$ (kV)', 'interpreter', 'latex')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

exportgraphics(fig, 'results/fig_case1_Vm.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Step 1: Detrending

% ---- Extracting voltage signals during the oscillation
y1 = Vm(1202:1802, 5);
N1 = size(y1, 1);
t1 = ts:ts:length(y1)*ts;

y2 = Vm(1202:1802,13);
N2 = size(y2, 1);
t2 = ts:ts:length(y2)*ts;

y1 = detrend(y1);
y2 = detrend(y2);

% ---- Plotting detrended signal during the event
width = 8;
height = 6;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')

subplot(211)
plot(t1, y1, 'color', '#191970')
xlim([0, max(t1)])
xlabel('Time after Event (s)')
ylabel('$\vert V \vert$ (kV)', 'interpreter', 'latex')
title('Voltage at Bus 5')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

subplot(212)
plot(t2, y2, 'color', '#191970')
xlim([0, max(t2)])
xlabel('Time after Event (s)')
ylabel('$\vert V \vert$ (kV)', 'interpreter', 'latex')
title('Voltage at Bus 13')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

exportgraphics(fig, 'results/fig_case1_Vm-event.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Step 2: Analyzing Signal Whiteness

% ---- Settings for autocorrelation plots
width = 8;
height = 6;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')

subplot(211)
autocorr(y1, 'NumLags', round(N1/4));
title('$\sigma_{y_1y_1}\left[\ell\right]$', 'interpreter', 'latex')

subplot(212)
autocorr(y2, 'NumLags', round(N2/4));
title('$\sigma_{y_2y_2}\left[\ell\right]$', 'interpreter', 'latex')

exportgraphics(fig, 'results/fig_case1_Vm-autocorr.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Step 3: Welch-based Power Spectral Density (PSD) estimation

% ---- Window-size
Nfft = [128, 256, 512];
linestyle = {'-', '-', '-', '-'};
color = {'#191970', '#6495ED', '#00BFFF', '#ADD8E6'};

% ---- Settings for PSD plots
width = 8;
height = 12;

fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
set(0, 'DefaultLineLineWidth', 2);
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultAxesFontName', 'Times')

for i=length(Nfft):-1:1

subplot(211)
% ---- PSD estimation using Welch method
[Pyy1, f] = pwelch(y1, Nfft(i), [], Nfft(i), fs);
plot(f, 10*log10(Pyy1), linestyle{i}, 'color', color{i}, 'DisplayName', ...
    "$N_{FFT}$ = " + string(Nfft(i)))
xlim([0, fs/2])
hold on

if i==1
    xlabel('$f$ [Hz]', 'interpreter', 'latex');
    ylabel('$10 \log_{10}\Vert \hat{S}_{y_1y_1}^2(f)\Vert$', ...
        'Interpreter', 'latex');
    title('Welch-based PSD Estimation of $y_1$ ($|V_{5}|$)', ...
        'Interpreter','latex');
    lgd = legend('interpreter', 'latex');
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
end

subplot(212)
[Pyy2, f] = pwelch(y2, Nfft(i), [], Nfft(i), fs);
plot(f, 10*log10(Pyy2), linestyle{i}, 'color', color{i}, 'DisplayName', ...
    "$N_{FFT}$ = " + string(Nfft(i)))
xlim([0, fs/2])
hold on

if i==1
    xlabel('$f$ [Hz]', 'interpreter', 'latex');
    ylabel('$10 \log_{10}\Vert \hat{S}_{y_2y_2}^2(f)\Vert$', ...
        'Interpreter', 'latex');
    title('Welch-based PSD Estimation of $y_2$ ($|V_{13}|$)', ...
        'Interpreter','latex');
    lgd = legend('interpreter', 'latex');
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
end

end

exportgraphics(fig, 'results/fig_case1_Vm-PSD.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

% -----------------------------------------------------------
% PSD for all voltage measurements
% -----------------------------------------------------------
clc
close all

% ---- Looping throughout all voltage measurements
for i=2:size(Vm, 2)

    % ---- Extracting and detrending signal
    y = Vm(1202:1802, i);
    N = size(y, 1);
    
    % ---- Skip any signal with non-numerical measurements
    if isnan(max(y))
        continue
    end
    
    % ---- Detrending
    y = detrend(y);
    
    % ---- Setting up plots
    width = 6;
    height = 6;

    fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
    set(0, 'DefaultLineLineWidth', 1.5);
    set(0,'DefaultAxesFontSize', 18)
    set(0,'DefaultTextFontSize', 18)
    set(0, 'DefaultTextFontName', 'Times')
    set(0, 'DefaultAxesFontName', 'Times')
    
    autocorr(y, 'NumLags', round(N/4));
    title_string = sprintf("{y_{%d} y_{%d}}", i, i);
    title(strcat('$\sigma_', title_string, '\left[\ell\right]$'), 'interpreter', 'latex')
    
    % ---- Exporting plot
    output_string = strcat('results/', ...
        sprintf('fig_case1_Vm_%d-autocorr.pdf', i));
    exportgraphics(fig, output_string, ...
        'ContentType', 'vector', 'BackGroundColor', 'none');

    % ---- PSD computation

    % ---- Window-size
    Nfft = [128, 256];
    linestyle = {'-', '-', '-', '-'};
    color = {'#191970', '#6495ED', '#00BFFF', '#ADD8E6'};
    
    % ---- PSD plot parameters
    width = 6;
    height = 6;

    fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
    set(0, 'DefaultLineLineWidth', 2);
    set(0,'DefaultAxesFontSize', 18)
    set(0,'DefaultTextFontSize', 18)
    set(0, 'DefaultTextFontName', 'Times')
    set(0, 'DefaultAxesFontName', 'Times')

    for j=length(Nfft):-1:1
        % ---- PSD estimation using Welch method
        [Pyy, f] = pwelch(y, Nfft(j), [], Nfft(j), fs);
        plot(f, 10*log10(Pyy), linestyle{j}, 'color', color{j}, 'DisplayName', ...
            "$N_{FFT}$ = " + string(Nfft(j)))
        xlim([0, fs/2])
        hold on
    
        if j==1

            title_string = sprintf("{y_{%d} y_{%d}}", i, i);

            xlabel('$f$ [Hz]', 'interpreter', 'latex');
            ylabel(strcat('$10 \log_{10}\Vert \hat{S}_', title_string, '^2(f)\Vert$'), ...
                'Interpreter', 'latex');

            title(sprintf('Welch-based PSD Estimation of $y_{%d}$ ($|V_{%d}|$)', i, i), ...
                'Interpreter','latex');
            lgd = legend('interpreter', 'latex');
            lgd.NumColumns = 4;
            lgd.Location = 'southoutside';
        end
    end

    % ---- Exporting plot
    output_string = strcat('results/', ...
        sprintf('fig_case1_Vm_%d-PSD.pdf', i));
    exportgraphics(fig, output_string, ...
        'ContentType', 'vector', 'BackGroundColor', 'none');

end