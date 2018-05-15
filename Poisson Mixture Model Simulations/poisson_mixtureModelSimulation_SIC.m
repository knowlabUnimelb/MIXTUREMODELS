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
hA = .05; % High Drift in channel A
hB = .05; % High Drift in channel B
lA = .02; % Low  Drift in channel A
lB = .02; % Low  Drift in channel B

pAB = 0;
pBA = pAB;

cA = 10; % Criterion A
cB = 10; % Criterion B

t = 0:15000; % Time vector to evaluate

pset = 0:.2:1;

drift = [hA, hB; hA lB; lA hB; lA lB];

%%
% for pidx = 1:size(pset,1)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Coactive cdfs
coactiveDrift.OR = drift - .0105;
coactiveDrift.AND = drift - .0105;
for j = 1:size(drift,1)
    coactive.cdf.OR(:,j) = fac_or([coactiveDrift.OR(j,:), 1, 1, cA, cB]', t)';
    coactive.pdf.OR(:,j) = diff([0; coactive.cdf.OR(:,j)]);
    
    coactive.cdf.AND(:,j) = fac_and([coactiveDrift.AND(j,:), 1, 1, cA, cB]', t)';
    coactive.pdf.AND(:,j) = diff([0; coactive.cdf.AND(:,j)]);
end
coactive.SIC.OR = sic(1-coactive.cdf.OR(:,[4 2 3 1]));
coactive.SIC.AND = sic(1-coactive.cdf.AND(:,[4 2 3 1]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Parallel cdfs
parallel.drift.AND = drift + .03;
parallel.drift.AND2 = drift - .01;
for j = 1:size(drift,1)
    parallel.cdf.OR(:,j) = fac_or([drift(j,:), 0, 0, cA, cB]', t)';
    parallel.pdf.OR(:,j) = diff([0; parallel.cdf.OR(:,j)]);
    
    parallel.cdf.AND(:,j) = fac_and([parallel.drift.AND(j,:), 0, 0, cA, cB]', t)';
    parallel.pdf.AND(:,j) = diff([0; parallel.cdf.AND(:,j)]);
    
    parallel.cdf.AND2(:,j) = fac_and([parallel.drift.AND2(j,:), 0, 0, cA, cB]', t)';
    parallel.pdf.AND2(:,j) = diff([0; parallel.cdf.AND2(:,j)]);
end

parallel.SIC.OR = sic(1-parallel.cdf.OR(:,[4 2 3 1]));
parallel.SIC.AND = sic(1-parallel.cdf.AND(:,[4 2 3 1]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get serial pdfs
pX = .5;
for j = 1:size(drift,1)
    A = singChan([drift(j,1), cA]', t)'; % PDF channel A
    B = singChan([drift(j,2), cB]', t)'; % PDF channel B
    
    pdfA = diff([0; A]);
    pdfB = diff([0; B]);
    
    serial.pdf.OR(:,j) = pX * pdfA + (1 - pX) * pdfB;
    serial.pdf.AND(:,j) = fftConv(pdfA, pdfB);
    
    serial.cdf.OR(:,j) =  cumsum(serial.pdf.OR(:,j)); % Self-terminating redundant target
    serial.cdf.AND(:,j) = cumsum(serial.pdf.AND(:,j));         % Exhaustive - convolve
end

serial.SIC.OR = sic(1-serial.cdf.OR(:,[4 2 3 1]));
serial.SIC.AND = sic(1-serial.cdf.AND(:,[4 2 3 1]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mix Models
cnt = 1;
for p = 0:.2:1;
    %%
    mixedCoaPar.pdf.OR(:,:,cnt) = p * coactive.pdf.OR + (1 - p) * parallel.pdf.OR;
    mixedCoaPar.cdf.OR(:,:,cnt) = cumsum(mixedCoaPar.pdf.OR(:,:,cnt));
    mixedCoaPar.S.OR(:,:,cnt)    = 1 - mixedCoaPar.cdf.OR(:,:,cnt);
    mixedCoaPar.S.OR(mixedCoaPar.S.OR <= 0) = 0;
    mixedCoaPar.SIC.OR(:,cnt) = sic(mixedCoaPar.S.OR(:,[4 2 3 1],cnt));
    
    %%
    mixedCoaPar.pdf.AND(:,:,cnt) = p * coactive.pdf.AND + (1 - p) * parallel.pdf.AND;
    mixedCoaPar.cdf.AND(:,:,cnt) = cumsum(mixedCoaPar.pdf.AND(:,:,cnt));
    mixedCoaPar.S.AND(:,:,cnt)    = 1 - mixedCoaPar.cdf.AND(:,:,cnt);
    mixedCoaPar.S.AND(mixedCoaPar.S.AND <= 0) = 0;
    mixedCoaPar.SIC.AND(:,cnt) = sic(mixedCoaPar.S.AND(:,[4 2 3 1],cnt));
    
    %%
    mixedSerPar.pdf.OR(:,:,cnt) = p * serial.pdf.OR + (1 - p) * parallel.pdf.OR;
    mixedSerPar.cdf.OR(:,:,cnt) = cumsum(mixedSerPar.pdf.OR(:,:,cnt));
    mixedSerPar.S.OR(:,:,cnt)    = 1 - mixedSerPar.cdf.OR(:,:,cnt);
    mixedCoaPar.S.OR(mixedCoaPar.S.OR <= 0) = 0;
    mixedSerPar.SIC.OR(:,cnt)   = sic(mixedSerPar.S.OR(:,[4 2 3 1], cnt));
    
    %%
    mixedSerPar.pdf.AND(:,:,cnt) = p * serial.pdf.AND + (1 - p) * parallel.pdf.AND2;
    mixedSerPar.cdf.AND(:,:,cnt) = cumsum(mixedSerPar.pdf.AND(:,:,cnt));
    mixedSerPar.S.AND(:,:,cnt)    = 1 - mixedSerPar.cdf.AND(:,:,cnt);
    mixedCoaPar.S.AND(mixedCoaPar.S.AND <= 0) = 0;
    mixedSerPar.SIC.AND(:,cnt)   = sic(mixedSerPar.S.AND(:,[4 2 3 1], cnt));
    
    cnt = cnt + 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
selection=1:10:3000;
colors = linspace(0, .8, 6)' * [1 1 1];
figure('WindowStyle', 'docked')
xlims = [0 1000];
ylims = [-.5 1];

subplot(2, 2, 1)
h = plot(selection', mixedCoaPar.SIC.OR(selection,:), '-');
set(gca, 'XLim', xlims, 'YLim', ylims, 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('SIC_{OR}(t)', 'FontSize', 12)
title('Parallel/Coactive', 'FontSize', 12)
legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

subplot(2, 2, 2)
h = plot(selection', mixedSerPar.SIC.OR(selection,:), '-');
set(gca, 'XLim', xlims, 'YLim', ylims, 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('SIC_{OR}(t)', 'FontSize', 12)
title('Parallel/Serial', 'FontSize', 12)
% legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

subplot(2, 2, 3)
h = plot(selection', mixedCoaPar.SIC.AND(selection,:), '-');
set(gca, 'XLim', xlims, 'YLim', ylims, 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('SIC_{AND}(t)', 'FontSize', 12)
title('Parallel/Coactive', 'FontSize', 12)
% legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

subplot(2, 2, 4)
h = plot(selection', mixedSerPar.SIC.AND(selection,:), '-');
set(gca, 'XLim', [0 2000], 'YLim', [-1 .5], 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('SIC_{AND}(t)', 'FontSize', 12)
title('Parallel/Serial', 'FontSize', 12)
% legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

toc;  