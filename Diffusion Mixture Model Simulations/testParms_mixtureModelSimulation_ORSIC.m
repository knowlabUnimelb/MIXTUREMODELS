% In v1, using the logical rule model set up, the model simulations did not
% necessary have the same mean RT or RT range resulting in some odd looking
% mixtures
% In v2, individual drift rates are tuned for each function to ensure that
% the mean RTs and variability are roughly commensurate
clear all
clc
close all

tic; 
settings          = struct();                % Set up GPU/CPU structure
settings.nGPUs    = 0;                       % Number of GPUs - if 0, use CPU
settings.gpuIds   = 0;                       % IDs for GPU
settings.nWalkers = 100000;                  % Number of GPU walkers
settings.seed     = 1;                       % Seed
[settings, ok] = checkSettings(settings,1); 

%% SFT Functions
sic    = @(S)(S(:,1) - S(:,2) - S(:,3) + S(:,4));        % LL, LH, HL, HH
capOR  = @(S)(log(S(:,1))./(log(S(:,2)) + log(S(:,3)))); % AB, A, B
capAND = @(F)((log(F(:,2)) + log(F(:,3)))./log(F(:,1))); % LL, A, B
ccf    = @(S)((log(S(:,1)) - log(S(:,2))) + (log(S(:,3)) - log(S(:,4)))); % AYh, AYl, XhB, XlB

% Set up possible parameters
a = [.5]; % Upper threshold
zr = .5;           % Where to start? .5 = halfway
sz = [.25]; % Width of uniform start point (as a proportion of a)
% v  - ?           % Drift rates
eta = [.1]; % vs 
% Ter = [.1 .2 .3];  % NDT
% sTer = [0 .05 .1 .2]; % NDT var

pset = allcomb(a, eta);
pset = [pset(:,1), pset(:,1) * sz(1), pset(:,2:end)];
nD = size(pset, 1);
for i = 1:numel(sz)-1
    pset = [pset; pset(1:nD, 1), pset(1:nD, 1) * sz(2), pset(1:nD, 3:end)]; 
end

%% Simulate AND SIC
% Drift rates were computed as a logit function of logical rules set up
coactive.driftrates = [8.16, 5.50, 5.50, 2.80]; % ABhh, ABhl, ABlh, AB
parallel.driftrates = [6.61, 1.67, 6.61, 1.67]; % H, L, H, L 
serial.driftrates   = parallel.driftrates; % H, L, H, L 
pX = .5; % Probability of sampling X first in a serial st model

coactive.stimnames     = {'ABhh', 'ABhl', 'ABlh', 'AB'};
parallel.dim.stimnames = {'Ah', 'A', 'Bh', 'B'};
serial.dim.stimnames   = parallel.dim.stimnames;
parallel.stimnames = coactive.stimnames;
serial.stimnames = parallel.stimnames;

%%
for pidx = 1:size(pset,1)
a = pset(pidx, 1); % Upper threshold
zr = .5;           % Where to start? .5 = halfway
sz = pset(pidx, 2); % Width of uniform start point (as a proportion of a)
eta = pset(pidx, 3); % vs 
Ter = 0; %pSet(pidx, 4);  % NDT
sTer = 0; %pSet(pidx, 5); % NDT var
    
%% Coactive LDMParms
coactive.diffusionparms          = struct();
coactive.diffusionparms.nStimuli = 4;
coactive.diffusionparms.a        = a; % (0 < a < 1)
coactive.diffusionparms.zr       = zr; % zr = z/a
coactive.diffusionparms.sz       = sz; % Width of uniform start point distribution U(zr_avg - sz/2, zr_avg + sz/2)
coactive.diffusionparms.v        = coactive.driftrates;
coactive.diffusionparms.eta      = eta * ones(1,4);
coactive.diffusionparms.Ter      = Ter;
coactive.diffusionparms.sTer     = sTer;

%% Parallel diffusionParms
parallel.diffusionparms          = struct();
parallel.diffusionparms.nStimuli = 4;
parallel.diffusionparms.a        = .5; % (0 < a < 1)
parallel.diffusionparms.zr       = zr; % zr = z/a
parallel.diffusionparms.sz       = .25; % Width of uniform start point distribution U(zr_avg - sz/2, zr_avg + sz/2) 
parallel.diffusionparms.v        = parallel.driftrates;
parallel.diffusionparms.eta      = .1 * ones(1,4);
parallel.diffusionparms.Ter      = Ter; % Mean Ter
parallel.diffusionparms.sTer     = sTer; % Width of uniform Ter

%% Serial diffusionParms
serial.diffusionparms          = struct();
serial.diffusionparms.nStimuli = 4;
serial.diffusionparms.a        = a; % (0 < a < 1)
serial.diffusionparms.zr       = zr; % zr = z/a
serial.diffusionparms.sz       = sz;
serial.diffusionparms.v        = serial.driftrates;
serial.diffusionparms.eta      = eta * ones(1,4);
serial.diffusionparms.Ter      = Ter; % Mean Ter
serial.diffusionparms.sTer     = sTer; % Width of uniform Ter

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Coactive diffusion pdfs
% tic;
[coactive.pdfs0, ok, extra, used] = LDMDist(coactive.diffusionparms, settings, 1);
coactive.diffusiont = toc;

% Process PDFs
coactive.pdfs.all = double(coactive.pdfs0)/double(settings.nWalkers);

coactive.pdfs.correct = coactive.pdfs.all(:,1:2:end); % LL, LH, HL, HH
coactive.pdfs.error   = coactive.pdfs.all(:,2:2:end); % LL, LH, HL, HH
coactive.cdfs.correct = cumsum(coactive.pdfs.correct,1);
coactive.S.correct = 1 - coactive.cdfs.correct;
coactive.orSIC = sic(coactive.S.correct(:, mstrfind(coactive.stimnames, {'AB', 'ABlh', 'ABhl', 'ABhh'})));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Parallel diffusion pdfs
%tic;
[parallel.dim.pdfs0, ok, extra, used] = LDMDist(parallel.diffusionparms, settings, 1);
parallel.diffusiont = toc;

% Process Parallel pdfs
parallel.dim.pdfs.all = double(parallel.dim.pdfs0)/double(settings.nWalkers);
parallel.dim.pdfs.correct = parallel.dim.pdfs.all(:,1:2:end); 
parallel.dim.pdfs.error   = parallel.dim.pdfs.all(:,2:2:end); 

parallel.dim.cdfs.correct = cumsum(parallel.dim.pdfs.all(:,1:2:end)); 
parallel.dim.cdfs.dim.error = cumsum(parallel.dim.pdfs.all(:,2:2:end)); 
parallel.dim.S.correct = 1 - parallel.dim.cdfs.correct; 
parallel.dim.S.error = 1 - parallel.dim.cdfs.correct; 

% Need to convert individual dimension pdfs to pdfs for each item
parallel.cdfs.correct = [1 - (parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'A'}))  .*  parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'B'}))) ,... % AB
                         1 - (parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'Ah'})) .*  parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'B'}))) ,... % ABhl
                         1 - (parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'A'}))  .*  parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'Bh'}))),... % ABLH
                         1 - (parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'Ah'})) .*  parallel.dim.S.correct(:, mstrfind(parallel.dim.stimnames, {'Bh'})))]; % LL
parallel.pdfs.correct = [zeros(1,4); diff(parallel.cdfs.correct)];
parallel.S.correct = 1 - parallel.cdfs.correct;
parallel.orSIC = sic(parallel.S.correct(:, mstrfind(parallel.stimnames, {'AB', 'ABlh', 'ABhl', 'ABhh'})));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get serial diffusion pdfs
%tic;
[serial.dim.pdfs0, ok, extra, used] = LDMDist(serial.diffusionparms, settings, 1);
% serial.pdfs0 = parallel.pdfs0;
% serial.diffusiont = toc;

% Process serial pdfs
serial.dim.pdfs.all = double(serial.dim.pdfs0)/double(settings.nWalkers);
serial.dim.pdfs.correct = serial.dim.pdfs.all(:,1:2:end); 
serial.dim.pdfs.error = serial.dim.pdfs.all(:,2:2:end); 

serial.dim.cdfs.correct = cumsum(serial.dim.pdfs.all(:,1:2:end)); 
serial.dim.cdfs.error = cumsum(serial.dim.pdfs.all(:,2:2:end)); 
serial.dim.S.correct = 1 - serial.dim.cdfs.correct; 
serial.dim.S.error = 1 - serial.dim.cdfs.correct; 

AB = pX * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'A'})) + (1 - pX) * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'B'}));
ABhl = pX * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'Ah'})) + (1 - pX) * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'B'}));
ABlh = pX * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'A'})) + (1 - pX) * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'Bh'}));
ABhh = pX * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'Ah'})) + (1 - pX) * serial.dim.pdfs.correct(:,mstrfind(serial.dim.stimnames, {'Bh'}));

% Need to convert individual dimension pdfs to pdfs for each item
serial.pdfs.correct = [AB, ABhl, ABlh, ABhh];
serial.cdfs.correct = cumsum(serial.pdfs.correct);
serial.S.correct = 1 - serial.cdfs.correct;
serial.orSIC = sic(serial.S.correct(:, mstrfind(coactive.stimnames, {'AB', 'ABlh', 'ABhl', 'ABhh'})));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mix Models
cnt = 1;
mixedCoaPar.stimnames = coactive.stimnames;
mixedSerPar.stimnames = coactive.stimnames;
mixedSerPar.dim.stimnames = parallel.dim.stimnames;
for p = 0:.2:1;
    %%
    mixedCoaPar.pdfs.correct(:,:,cnt) = p * coactive.pdfs.correct + (1 - p) * parallel.pdfs.correct;
    mixedCoaPar.cdfs.correct(:,:,cnt) = cumsum(mixedCoaPar.pdfs.correct(:,:,cnt));
    mixedCoaPar.S.correct(:,:,cnt)    = 1 - mixedCoaPar.cdfs.correct(:,:,cnt); 
    mixedCoaPar.S.correct(mixedCoaPar.S.correct < 0) = 0;
    
    mixedCoaPar.orSIC(:,cnt) = sic(mixedCoaPar.S.correct(:,mstrfind(mixedCoaPar.stimnames, {'AB', 'ABlh', 'ABhl', 'ABhh'}), cnt));

    %% 
    mixedSerPar.dim.pdfs.correct(:,:,cnt) = p * serial.dim.pdfs.correct + (1 - p) * parallel.dim.pdfs.correct;
    mixedSerPar.dim.cdfs.correct(:,:,cnt) = cumsum(mixedSerPar.dim.pdfs.correct(:,:,cnt));
    mixedSerPar.dim.S.correct(:,:,cnt)    = 1 - mixedSerPar.dim.cdfs.correct(:,:,cnt);    
    
    mixedSerPar.pdfs.correct(:,:,cnt) = p * serial.pdfs.correct + (1 - p) * parallel.pdfs.correct;
    mixedSerPar.cdfs.correct(:,:,cnt) = cumsum(mixedSerPar.pdfs.correct(:,:,cnt));
    mixedSerPar.S.correct(:,:,cnt)    = 1 - mixedSerPar.cdfs.correct(:,:,cnt);

    mixedSerPar.orSIC(:,cnt)   = sic(mixedSerPar.S.correct(:,mstrfind(mixedSerPar.stimnames, {'AB', 'ABlh', 'ABhl', 'ABhh'}), cnt));

    cnt = cnt + 1; 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 
selection=1:10:3000;
%%
figure('WindowStyle', 'docked')
fieldnames = {'orSIC'};
xlims = [0 250];
ylims = [-.5 1];
colors = linspace(0, .8, 6)' * [1 1 1];
for i = 1%5
    subplot(1, 2, i)
   h = plot(selection', mixedCoaPar.(fieldnames{i})(selection,:), '-');
       for j = 1:size(colors, 1)
        set(h(j), 'Color', colors(j,:), 'LineWidth', 2)
    end
    set(gca, 'XLim', xlims(i,:), 'YLim', ylims(i,:), 'FontSize', 10)
    xt = get(gca,'XTick');
    set(gca, 'XTickLabel', num2str(xt' * 10))
    xlabel('t (msec)', 'FontSize', 12)
    ylabel('SIC_{OR}(t)', 'FontSize', 12)
    title('Parallel/Coactive', 'FontSize', 12)
    legend(num2str((0:.2:1)'))
    
    subplot(1, 2, i+1)
   h =  plot(selection', mixedSerPar.(fieldnames{i})(selection,:), '-');
       for j = 1:size(colors, 1)
        set(h(j), 'Color', colors(j,:), 'LineWidth', 2)
    end
    set(gca, 'XLim', xlims(i,:), 'YLim', ylims(i,:), 'FontSize', 10)
    xt = get(gca,'XTick');
    set(gca, 'XTickLabel', num2str(xt' * 10))
    xlabel('t (msec)', 'FontSize', 12)
    ylabel('SIC_{OR}(t)', 'FontSize', 12)
    title('Parallel/Serial', 'FontSize', 12)
    legend(num2str((0:.2:1)'))
end
% % displayTime(toc)
toc
end