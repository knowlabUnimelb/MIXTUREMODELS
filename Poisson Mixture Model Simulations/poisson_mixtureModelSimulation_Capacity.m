clear all
clc
close all

tic

%% SFT Functions
sic    = @(S)(S(:,1) - S(:,2) - S(:,3) + S(:,4));        % LL, LH, HL, HH
capOR  = @(S)(log(S(:,1))./(log(S(:,2)) + log(S(:,3)))); % AB, A, B
capAND = @(F)((log(F(:,2)) + log(F(:,3)))./log(F(:,1))); % LL, A, B
ccf    = @(S)((log(S(:,1)) - log(S(:,2))) + (log(S(:,3)) - log(S(:,4)))); % AYh, AYl, XhB, XlB

% Set up possible parameters
A = .05; % High Drift in channel A
B = .05; % High Drift in channel B
AB = [A B];

pAB = 0;
pBA = pAB;

cA = 10; % Criterion A
cB = 10; % Criterion B

t = 0:15000; % Time vector to evaluate

pset = 0:.2:1;


%%
% for pidx = 1:size(pset,1)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Coactive cdfs
coactive.cdf.OR(:,1) = fac_or([AB, 1, 1, cA, cB]', t)'; % AB
coactive.pdf.OR(:,1) = diff([0; coactive.cdf.OR(:,1)]);
      
coactive.cdf.OR(:,2) = singChan([A, cA]', t)';
coactive.pdf.OR(:,2) = diff([0; coactive.cdf.OR(:,2)]);

coactive.cdf.OR(:,3) = singChan([B, cB]', t)';
coactive.pdf.OR(:,3) = diff([0; coactive.cdf.OR(:,2)]);
    
coactive.cdf.AND(:,1) = fac_and([AB, 1, 1, cA, cB]', t)'; % AB
coactive.pdf.AND(:,1) = diff([0; coactive.cdf.AND(:,1)]);

coactive.cdf.AND(:,2) = singChan([A, cA]', t)';
coactive.pdf.AND(:,2) = diff([0; coactive.cdf.OR(:,2)]);

coactive.cdf.AND(:,3) = singChan([B, cB]', t)';
coactive.pdf.AND(:,3) = diff([0; coactive.cdf.OR(:,2)]);

coactive.cap.OR = capOR(1-coactive.cdf.OR(:,[1 2 3]));
coactive.cap.AND = capAND(coactive.cdf.AND(:,[1 2 3]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Parallel cdfs
parallel.cdf.OR(:,1) = fac_or([AB, 0, 0, cA, cB]', t)'; % AB
parallel.pdf.OR(:,1) = diff([0; parallel.cdf.OR(:,1)]);
      
parallel.cdf.OR(:,2) = singChan([A, cA]', t)';
parallel.pdf.OR(:,2) = diff([0; parallel.cdf.OR(:,2)]);

parallel.cdf.OR(:,3) = singChan([B, cB]', t)';
parallel.pdf.OR(:,3) = diff([0; parallel.cdf.OR(:,2)]);
    
parallel.cdf.AND(:,1) = fac_and([AB, 0, 0, cA, cB]', t)';  % AB
parallel.pdf.AND(:,1) = diff([0; parallel.cdf.AND(:,1)]);

parallel.cdf.AND(:,2) = singChan([A, cA]', t)';
parallel.pdf.AND(:,2) = diff([0; parallel.cdf.OR(:,2)]);

parallel.cdf.AND(:,3) = singChan([B, cB]', t)';
parallel.pdf.AND(:,3) = diff([0; parallel.cdf.OR(:,2)]);

parallel.cap.OR = capOR(1-parallel.cdf.OR(:,[1 2 3]));
parallel.cap.AND = capAND(parallel.cdf.AND(:,[1 2 3]));

parallel.cap.OR = capOR(1-parallel.cdf.OR(:,[1 2 3]));
parallel.cap.AND = capAND(parallel.cdf.AND(:,[1 2 3]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get serial pdfs
pX = .5;

A = singChan([AB(1), cA]', t)'; % PDF channel A
B = singChan([AB(2), cB]', t)'; % PDF channel B

pdfA = diff([0; A]);
pdfB = diff([0; B]);
    
serial.pdf.OR = pX * pdfA + (1 - pX) * pdfB;
serial.pdf.AND = fftConv(pdfA, pdfB);
 
serial.cdf.OR  =  cumsum(serial.pdf.OR); % Self-terminating redundant target
serial.cdf.AND = cumsum(serial.pdf.AND); % Exhaustive - convolve

serial.cdf.OR(:,2) = A;
serial.cdf.OR(:,3) = B;
serial.cdf.AND(:,2) = A;
serial.cdf.AND(:,3) = B;
serial.pdf.OR(:,2) = pdfA;
serial.pdf.OR(:,3) = pdfB;
serial.pdf.AND(:,2) = pdfA;
serial.pdf.AND(:,3) = pdfB;

parallel.cap.OR = capOR(1-serial.cdf.OR(:,[1 2 3]));
parallel.cap.AND = capAND(serial.cdf.AND(:,[1 2 3]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mix Models
cnt = 1;
for p = 0:.2:1;
    %%
    mixedCoaPar.pdf.OR(:,:,cnt) = p * coactive.pdf.OR + (1 - p) * parallel.pdf.OR;
    mixedCoaPar.cdf.OR(:,:,cnt) = cumsum(mixedCoaPar.pdf.OR(:,:,cnt));
    mixedCoaPar.S.OR(:,:,cnt)    = 1 - mixedCoaPar.cdf.OR(:,:,cnt);
    mixedCoaPar.S.OR(mixedCoaPar.S.OR <= 0) = 0;
    mixedCoaPar.cap.OR(:,cnt) = capOR(mixedCoaPar.S.OR(:,[1 2 3],cnt));
    
    %%
    mixedCoaPar.pdf.AND(:,:,cnt) = p * coactive.pdf.AND + (1 - p) * parallel.pdf.AND;
    mixedCoaPar.cdf.AND(:,:,cnt) = cumsum(mixedCoaPar.pdf.AND(:,:,cnt));
    mixedCoaPar.S.AND(:,:,cnt)    = 1 - mixedCoaPar.cdf.AND(:,:,cnt);
    mixedCoaPar.S.AND(mixedCoaPar.S.AND <= 0) = 0;
    mixedCoaPar.cap.AND(:,cnt) = capAND(mixedCoaPar.cdf.AND(:,[1 2 3],cnt));
    
    %%
    mixedSerPar.pdf.OR(:,:,cnt) = p * serial.pdf.OR + (1 - p) * parallel.pdf.OR;
    mixedSerPar.cdf.OR(:,:,cnt) = cumsum(mixedSerPar.pdf.OR(:,:,cnt));
    mixedSerPar.S.OR(:,:,cnt)    = 1 - mixedSerPar.cdf.OR(:,:,cnt);
    mixedCoaPar.S.OR(mixedCoaPar.S.OR <= 0) = 0;
    mixedSerPar.cap.OR(:,cnt)   = capOR(mixedSerPar.S.OR(:,[1 2 3], cnt));
    
    %%
    mixedSerPar.pdf.AND(:,:,cnt) = p * serial.pdf.AND + (1 - p) * parallel.pdf.AND;
    mixedSerPar.cdf.AND(:,:,cnt) = cumsum(mixedSerPar.pdf.AND(:,:,cnt));
    mixedSerPar.S.AND(:,:,cnt)    = 1 - mixedSerPar.cdf.AND(:,:,cnt);
    mixedCoaPar.S.AND(mixedCoaPar.S.AND <= 0) = 0;
    mixedSerPar.cap.AND(:,cnt)   = capAND(mixedSerPar.cdf.AND(:,[1 2 3], cnt));
    
    cnt = cnt + 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
selection=1:10:3000;
colors = linspace(0, .8, 6)' * [1 1 1];

figure('WindowStyle', 'docked')
xlims = [0 500];
ylims = [0 1.1];

subplot(2, 2, 1)
h = plot(selection', mixedCoaPar.cap.OR(selection,:), '-');
set(gca, 'XLim', xlims, 'YLim', [0 15], 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('C_{OR}(t)', 'FontSize', 12)
title('Parallel/Coactive', 'FontSize', 12)
legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

subplot(2, 2, 2)
h = plot(selection', mixedSerPar.cap.OR(selection,:), '-');
set(gca, 'XLim', xlims, 'YLim', ylims, 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('C_{OR}(t)', 'FontSize', 12)
title('Parallel/Serial', 'FontSize', 12)
legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

subplot(2, 2, 3)
h = plot(selection', mixedCoaPar.cap.AND(selection,:), '-');
set(gca, 'XLim', xlims, 'YLim', [0 15], 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('C_{AND}(t)', 'FontSize', 12)
title('Parallel/Coactive', 'FontSize', 12)
legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

subplot(2, 2, 4)
h = plot(selection', mixedSerPar.cap.AND(selection,:), '-');
set(gca, 'XLim', xlims, 'YLim', ylims, 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('C_{AND}(t)', 'FontSize', 12)
title('Parallel/Serial', 'FontSize', 12)
legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

toc;  