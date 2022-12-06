clc
clear
close all

% Loading given data
load('data\sysid_power_case2.mat');

%% Visualization

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

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case2_Vm-time-plot.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Data Preprocessing

y1 = Vm(1:end-800, 7);
N1 = size(y1, 1);
t1 = ts:ts:length(y1)*ts;

y2 = Vm(1:end-800, 10);
N2 = size(y2, 1);
t2 = ts:ts:length(y2)*ts;

y1 = detrend(y1);
y2 = detrend(y2);

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
title('Voltage at Bus 7')
xlabel('Time (s)')
ylabel('$|V|$ (kV)', 'interpreter', 'latex')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

subplot(212)
plot(t2, y2,  'color', '#191970')
xlim([0, max(t2)])
title('Voltage at Bus 10')
xlabel('Time (s)')
ylabel('$|V|$ (kV)', 'interpreter', 'latex')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case2_Vm-detrend2.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Autocorrelation Plot

% ================================
% Analyzing whiteness of signals via autocorrelation
% ================================
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

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case2_Vm-autocorr.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Welch-based PSD

% ================================
% Plotting Welch-based PSD of output
% ================================

% Window-size
Nfft = [1024, 2048, 5012];
linestyle = {'-', '-', '-', '-'};
color = {'#191970', '#6495ED', '#00BFFF', '#ADD8E6'};

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
% PSD estimation using Welch method
[Pyy1, f] = pwelch(y1, Nfft(i), [], Nfft(i), fs);
plot(f, 10*log10(Pyy1), linestyle{i}, 'color', color{i}, 'DisplayName', ...
    "$N_{FFT}$ = " + string(Nfft(i)))
xlim([0, fs/2])
hold on

if i==1
    xlabel('$f$ [Hz]', 'interpreter', 'latex');
    ylabel('$10 \log_{10}\Vert \hat{S}_{y_1y_1}^2(f)\Vert$', ...
        'Interpreter', 'latex');
    title('Welch-based PSD Estimation of $y_1$ ($|V_{7}|$)', ...
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
    title('Welch-based PSD Estimation of $y_2$ ($|V_{10}|$)', ...
        'Interpreter','latex');
    lgd = legend('interpreter', 'latex');
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
end

end

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case2_Vm-PSD.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');