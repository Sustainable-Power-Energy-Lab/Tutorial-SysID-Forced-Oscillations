clc
clear
close all

% Loading given data
load('data\sysid_power_case1.mat');

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
plot(time, Freq(:, 5), 'color', '#191970')
xlim([0, max(time)])
title('Frequency at Bus 5')
xlabel('Time (s)')
ylabel('$f$ (Hz)', 'interpreter', 'latex')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

subplot(212)
plot(time, Freq(:, 13),  'color', '#191970')
xlim([0, max(time)])
title('Frequency at Bus 13')
xlabel('Time (s)')
ylabel('$f$ (Hz)', 'interpreter', 'latex')
xline(t_event, 'color', 'red', 'linewidth', 2, 'linestyle', '--')

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case1_f.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Analyzing Waveform During Event

y1 = Freq(1202:1802, 5);
N1 = size(y1, 1);
t1 = ts:ts:length(y1)*ts;

y2 = Freq(1202:1802,13);
N2 = size(y2, 1);
t2 = ts:ts:length(y2)*ts;

y1 = detrend(y1);
y2 = detrend(y2);

% ===============================
% Time-domain plot of the signals (during the event)
% ===============================

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
ylabel('$f$ (Hz)', 'interpreter', 'latex')
title('Frequency at Bus 5')

subplot(212)
plot(t2, y2, 'color', '#191970')
xlim([0, max(t2)])
xlabel('Time after Event (s)')
ylabel('$f$ (Hz)', 'interpreter', 'latex')
title('Frequency at Bus 5')

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case1_f-event.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Autocorrelation Plot

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

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case1_f-autocorr.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Welch-based PSD

% ================================
% Plotting Welch-based PSD of input and output
% ================================

% Window-size
Nfft = [128, 256, 512];
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
    title('Welch-based PSD Estimation of $y_1$ ($f_{5}$)', ...
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
    title('Welch-based PSD Estimation of $y_2$ ($f_{13}$)', ...
        'Interpreter','latex');
    lgd = legend('interpreter', 'latex');
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
end

end

exportgraphics(fig, 'C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\fig_case1_f-PSD.pdf', ...
    'ContentType', 'vector', 'BackGroundColor', 'none');

%% Analysis for all Measurements

clc
close all

% Looping throughout all voltage measurements
for i=2:size(Freq, 2)

    % ===============================
    % Preprocessing signals
    % ===============================
    y = Freq(1202:1802, i);
    N = size(y, 1);

    if isnan(max(y))
        continue
    end

    y = detrend(y);

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
    
    % Exporting plot
    output_string = strcat('C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\', ...
        sprintf('fig_case1_f_%d-autocorr.pdf', i));
%     exportgraphics(fig, output_string, ...
%         'ContentType', 'vector', 'BackGroundColor', 'none');

    % =============
    % Welch-based PSD
    % =============

    % Window-size
    Nfft = [128, 256];
    linestyle = {'-', '-', '-', '-'};
    color = {'#191970', '#6495ED', '#00BFFF', '#ADD8E6'};
    
    width = 6;
    height = 6;

    fig = figure('Renderer', 'opengl', 'Position', [10 10 width*100 height*100]);
    set(0, 'DefaultLineLineWidth', 2);
    set(0,'DefaultAxesFontSize', 18)
    set(0,'DefaultTextFontSize', 18)
    set(0, 'DefaultTextFontName', 'Times')
    set(0, 'DefaultAxesFontName', 'Times')

    for j=length(Nfft):-1:1
        % PSD estimation using Welch method
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

            title(sprintf('Welch-based PSD Estimation of $y_{%d}$ ($|I_{%d}|$)', i, i), ...
                'Interpreter','latex');
            lgd = legend('interpreter', 'latex');
            lgd.NumColumns = 4;
            lgd.Location = 'southoutside';
        end
    end

    % Exporting plot
    output_string = strcat('C:\Users\Sergio\Insync\sergio.dorado.rojas@gmail.com\Dropbox\Apps\Overleaf\sysid_project_presentation\figs\', ...
        sprintf('fig_case1_Im_%d-PSD.pdf', i));
%     exportgraphics(fig, output_string, ...
%         'ContentType', 'vector', 'BackGroundColor', 'none');

end