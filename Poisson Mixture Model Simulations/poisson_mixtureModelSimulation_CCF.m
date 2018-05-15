% Poisson channel mixture model simulation
% For "Nice Guys Bat for Both Teams" paper
% 170228 DL
clear all
clc
close all

%% SFT Functions
sic    = @(S)(S(:,1) - S(:,2) - S(:,3) + S(:,4));        % LL, LH, HL, HH
capOR  = @(S)(log(S(:,1))./(log(S(:,2)) + log(S(:,3)))); % AB, A, B
capAND = @(F)((log(F(:,2)) + log(F(:,3)))./log(F(:,1))); % LL, A, B
ccf    = @(S)((log(S(:,1)) - log(S(:,2))) + (log(S(:,3)) - log(S(:,4)))); % AYh, AYl, XhB, XlB

%% Parameters
% Channel "drift"
vA  = .02;  % Target Channel A (target channels trigger termination of the process)
vB  = .02;  % Target Channel B (target channels trigger termination of the process)
vYL = .02; %.015;  % Conflict channel YL (low salience conflict)
vYH = .05; %.025; % Conflict channel YH (high salience conflict)
vXL = .02; %.015;  % Conflict channel XL (low salience conflict)
vXH = .05; %.025; % Conflict channel XH (high salience conflict)
drift = [vA vYH; vA vYL; vB vXH; vB vXL];

% Thresholds
cA = 3; % Criterion A
cB = 3; % Criterion B

%%
t = 0:15000; % Time vector to evaluate

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Coactive cdfs
for j = 1:size(drift,1)
%     coactive.cdf.OR(:,j) = fac_or([drift(j,:), 1, 1, cA, cB]', t)';  % Use a facilitatory parallel conflict model with complete cross-talk
    coactive.cdf.OR(:,j) = inh_conf([drift(j,:), 1, 1, cA, cB]', t)';
    coactive.pdf.OR(:,j) = diff([0; coactive.cdf.OR(:,j)]);            % differentiate cdf to get pdf
    
%     coactive.cdf.AND(:,j) = fac_and([drift(j,:), 1, 1, cA, cB]', t)'; % Use a facilitatory parallel conflict model with complete cross-talk
    coactive.cdf.AND(:,j) = inh_conf([drift(j,:), 1, 1, cA, cB]', t)';
    coactive.pdf.AND(:,j) = diff([0; coactive.cdf.AND(:,j)]);          % differentiate cdf to get pdf
end
coactive.CCF.OR  = ccf(1-coactive.cdf.OR(:,[1 2 3 4]));
coactive.CCF.AND = ccf(1-coactive.cdf.AND(:,[1 2 3 4]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Parallel cdfs
for j = 1:size(drift,1)
    parallel.cdf.OR(:,j) = inh_conf([drift(j,:), 0, 0, cA, cB]', t)'; % Use a facilitatory parallel model with no cross-talk
    parallel.pdf.OR(:,j) = diff([0; parallel.cdf.OR(:,j)]);
    
    parallel.cdf.AND(:,j) = inh_conf([drift(j,:), 0, 0, cA, cB]', t)'; % No cross-talk so just use fac_conf
    parallel.pdf.AND(:,j) = diff([0; parallel.cdf.AND(:,j)]);
end

parallel.CCF.OR = sic(1-parallel.cdf.OR(:,[1 2 3 4]));
parallel.CCF.AND = sic(1-parallel.cdf.AND(:,[1 2 3 4]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get serial pdfs
pX = .5;

A = singChan([vA, cA]', t)'; % PDF channel A
B = singChan([vB, cB]', t)'; % PDF channel B
YL = singChan([vYL, cA]', t)'; % PDF channel YL
YH = singChan([vYH, cB]', t)'; % PDF channel YH
XL = singChan([vXL, cA]', t)'; % PDF channel XL
XH = singChan([vXH, cB]', t)'; % PDF channel XH

pdfA = diff([0; A]);
pdfB = diff([0; B]);
pdfYL = diff([0; YL]);
pdfYH = diff([0; YH]);
pdfXL = diff([0; XL]);
pdfXH = diff([0; XH]);

AYh = pX * pdfA + (1 - pX) * fftConv(pdfA, pdfYH);
AYl = pX * pdfA + (1 - pX) * fftConv(pdfA, pdfYL);
XhB = pX * pdfB + (1 - pX) * fftConv(pdfB, pdfXH);
XlB = pX * pdfB + (1 - pX) * fftConv(pdfB, pdfXL);

exAYh = fftConv(pdfA, pdfYH);
exAYl = fftConv(pdfA, pdfYL);
exXhB = fftConv(pdfB, pdfXH);
exXlB = fftConv(pdfB, pdfXL);

serial.pdf.OR = [AYh, AYl, XhB, XlB];
serial.pdf.AND = [exAYh, exAYl, exXhB, exXlB];

for j = 1:4
serial.cdf.OR(:,j) =  cumsum(serial.pdf.OR(:,j)); % Self-terminating redundant target
serial.cdf.AND(:,j) = cumsum(serial.pdf.AND(:,j));         % Exhaustive - convolve
end
serial.CCF.OR = ccf(1-serial.cdf.OR(:,[1 2 3 4]));
serial.CCF.AND = ccf(1-serial.cdf.AND(:,[1 2 3 4]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mix Models
cnt = 1;
for p = 0:.2:1;
    %%
    mixedCoaPar.pdf.OR(:,:,cnt) = p * coactive.pdf.OR + (1 - p) * parallel.pdf.OR;
    mixedCoaPar.cdf.OR(:,:,cnt) = cumsum(mixedCoaPar.pdf.OR(:,:,cnt));
    mixedCoaPar.S.OR(:,:,cnt)    = 1 - mixedCoaPar.cdf.OR(:,:,cnt);
    mixedCoaPar.S.OR(mixedCoaPar.S.OR <= 0) = 0;
    mixedCoaPar.CCF.OR(:,cnt) = ccf(mixedCoaPar.S.OR(:,[1 2 3 4],cnt));
    
    %%
    mixedCoaPar.pdf.AND(:,:,cnt) = p * coactive.pdf.AND + (1 - p) * parallel.pdf.AND;
    mixedCoaPar.cdf.AND(:,:,cnt) = cumsum(mixedCoaPar.pdf.AND(:,:,cnt));
    mixedCoaPar.S.AND(:,:,cnt)    = 1 - mixedCoaPar.cdf.AND(:,:,cnt);
    mixedCoaPar.S.AND(mixedCoaPar.S.AND <= 0) = 0;
    mixedCoaPar.CCF.AND(:,cnt) = ccf(mixedCoaPar.S.AND(:,[1 2 3 4],cnt));
    
    %%
    mixedSerPar.pdf.OR(:,:,cnt) = p * serial.pdf.OR + (1 - p) * parallel.pdf.OR;
    mixedSerPar.cdf.OR(:,:,cnt) = cumsum(mixedSerPar.pdf.OR(:,:,cnt));
    mixedSerPar.S.OR(:,:,cnt)    = 1 - mixedSerPar.cdf.OR(:,:,cnt);
    mixedCoaPar.S.OR(mixedCoaPar.S.OR <= 0) = 0;
    mixedSerPar.CCF.OR(:,cnt)   = ccf(mixedSerPar.S.OR(:,[1 2 3 4], cnt));
    
    %%
    mixedSerPar.pdf.AND(:,:,cnt) = p * serial.pdf.AND + (1 - p) * parallel.pdf.AND;
    mixedSerPar.cdf.AND(:,:,cnt) = cumsum(mixedSerPar.pdf.AND(:,:,cnt));
    mixedSerPar.S.AND(:,:,cnt)    = 1 - mixedSerPar.cdf.AND(:,:,cnt);
    mixedCoaPar.S.AND(mixedCoaPar.S.AND <= 0) = 0;
    mixedSerPar.CCF.AND(:,cnt)   = ccf(mixedSerPar.S.AND(:,[1 2 3 4], cnt));
    
    cnt = cnt + 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
selection=1:10:3000;
colors = linspace(0, .8, 6)' * [1 1 1];
figure('WindowStyle', 'docked')
xlims = [0 1000];
ylims = [-.5 1];

subplot(1,2, 1)
h=plot(selection', mixedCoaPar.CCF.OR(selection,:), '-');
set(gca, 'XLim', [0 1500], 'YLim', [-.1 2], 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('CCF(t)', 'FontSize', 12)
title('Parallel/Coactive', 'FontSize', 12)
legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

subplot(1,2, 2)
h=plot(selection', mixedSerPar.CCF.OR(selection,:), '-');
set(gca, 'XLim', [0 1500], 'YLim', [-2 .1], 'FontSize', 10)
xlabel('t (msec)', 'FontSize', 12)
ylabel('CCF(t)', 'FontSize', 12)
title('Parallel/Serial', 'FontSize', 12)
% legend(num2str((0:.2:1)'))
for i = 1:size(colors, 1)
    set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
end

% subplot(2, 2, 3)
% h=plot(selection', mixedCoaPar.CCF.AND(selection,:), '-');
% set(gca, 'XLim', xlims, 'YLim', ylims, 'FontSize', 10)
% xlabel('t (msec)', 'FontSize', 12)
% ylabel('CCF(t)', 'FontSize', 12)
% title('Coactive/Parallel Exhaustive', 'FontSize', 12)
% % legend(num2str((0:.2:1)'))
% for i = 1:size(colors, 1)
%     set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
% end

% subplot(2, 2, 4)
% h=plot(selection', mixedSerPar.CCF.AND(selection,:), '-');
% set(gca, 'XLim', [0 2000], 'YLim', [-1 .5], 'FontSize', 10)
% xlabel('t (msec)', 'FontSize', 12)
% ylabel('CCF(t)', 'FontSize', 12)
% title('Serial/Parallel Exhaustive', 'FontSize', 12)
% % legend(num2str((0:.2:1)'))
% for i = 1:size(colors, 1)
%     set(h(i), 'Color', colors(i,:), 'LineWidth', 2)
% end

% toc;  