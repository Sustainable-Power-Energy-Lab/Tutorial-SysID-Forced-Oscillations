clc
close all

% Time vector
t = ScopeData1.time;

% Signal values
x1 = ScopeData1.signals(1).values; % current
x2 = ScopeData1.signals(2).values; % voltage
x3 = ScopeData1.signals(3).values; % p and q

ia = x2(:, 1);
ib = x2(:, 2);
ic = x2(:, 3);

va = x1(:, 1);
vb = x1(:, 2);
vc = x1(:, 3);

P = x3(:, 1);
Q = x3(:, 2);

header = {'t', 'va', 'vb', 'vc', 'ia', 'ib', 'ic', 'P', 'Q'};
header= strjoin(header, ',');
fid = fopen('Converter_VIPQ-Data.csv','w'); 
fprintf(fid,'%s\n',header);
fclose(fid);
rol
dlmwrite('Converter_VIPQ-Data.csv', [t va vb vc ia ib ic P Q], '-append')