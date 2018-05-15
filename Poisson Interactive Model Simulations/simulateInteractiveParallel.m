%% Simulate Interactive Parallel Model and plot results
% Top row are the results for the facilitatory parallel model
% Bottom row are the results for the inhibitory parallel model
%
% Columns should be: SIC_OR, SIC_AND, C_OR, C_AND, CCF

clear all
clc
close all
tic

runSICf = true;
runCf = true;
runCCFf = true;
runSICi = true;
runCi = true;
runCCFi = true;

colors = linspace(0, .6, 5)' * [1 1 1];

%%
hA = .05; % High Drift in channel A
hB = .05; % High Drift in channel B
lA = .02; % Low  Drift in channel A
lB = .02; % Low  Drift in channel B

pAB = [0 .25 .5 .75 1];
pBA = pAB;

cA = 10; % Criterion A
cB = 10; % Criterion B

t = 0:15000; % Time vector to evaluate

drift = [hA, hB; hA lB; lA hB; lA lB];

%% Get cdfs and SICs
if runSICf
    for i = 1:numel(pAB)
        for j = 1:size(drift,1)
            facilitatory.cdf.OR(:,j,i) = fac_or([drift(j,:), pAB(i), pBA(i), cA, cB]', t)';
            facilitatory.cdf.AND(:,j,i) = fac_and([drift(j,:), pAB(i), pBA(i), cA, cB]', t)';
        end
        facilitatory.SIC.OR(:,i) =...
            (1-facilitatory.cdf.OR(:,4,i)) - (1-facilitatory.cdf.OR(:,2,i)) -...
            (1-facilitatory.cdf.OR(:,3,i)) + (1-facilitatory.cdf.OR(:,1,i));
        facilitatory.SIC.AND(:,i) =...
            (1-facilitatory.cdf.AND(:,4,i)) - (1-facilitatory.cdf.AND(:,2,i)) -...
            (1-facilitatory.cdf.AND(:,3,i)) + (1-facilitatory.cdf.AND(:,1,i));
    end
    
    %% Plot SIC-OR results
    subplot(2,5,1)
    line([min(t), max(t)], [0 0], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t, facilitatory.SIC.OR);
    legend(h, num2str(pAB'))
    set(gca, 'XLim', [0 1000], 'YLim', [-.75 1.25], 'FontSize', 12, 'XTick', [])
    ylabel('SIC_{OR}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
    
    %% Plot SIC-AND results
    subplot(2,5,2)
    line([min(t), max(t)], [0 0], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t, facilitatory.SIC.AND);
    set(gca, 'XLim', [0 1000], 'YLim', [-.75 1.25], 'FontSize', 12, 'XTick', [])
    ylabel('SIC_{AND}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    toc
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
end

%% Get Capacity
A = .02;
B = .02;
if runCf
    for i = 1:numel(pAB)
        facilitatory.ccdf.OR(:,1,i) = fac_or([A, B, pAB(i), pBA(i), cA, cB]', t)';
        facilitatory.ccdf.OR(:,2,i) = ind_sing([A, cA]', t)';
        facilitatory.ccdf.OR(:,3,i) = ind_sing([B, cB]', t)';
        
        facilitatory.ccdf.AND(:,1,i) = fac_and([A, B, pAB(i), pBA(i), cA, cB]', t)';
        facilitatory.ccdf.AND(:,2,i) = facilitatory.ccdf.OR(:,2,i);
        facilitatory.ccdf.AND(:,3,i) = facilitatory.ccdf.OR(:,3,i);
        
        facilitatory.CAP.OR(:,i) =...
            log(1 - facilitatory.ccdf.OR(:,1,i))./...
            log((1 - facilitatory.ccdf.OR(:,2,i)) .* (1 - facilitatory.ccdf.OR(:,3,i)));      % AB/(A + B);
        facilitatory.CAP.AND(:,i) =...
            log((facilitatory.ccdf.AND(:,2,i)) .* (facilitatory.ccdf.AND(:,3,i)))./...
            log(facilitatory.ccdf.AND(:,1,i));      % (A + B)/AB;
    end
    
    %% Plot Capacity OR results
    subplot(2,5,3)
    line([min(t), max(t)], [1 1], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t,facilitatory.CAP.OR);
    set(gca, 'XLim', [0 1500], 'YLim', [-.5 6], 'FontSize', 12, 'XTick', [])
    ylabel('C_{OR}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    % supertitle('OR Capacity')
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
    
    %% Plot Capacity AND results
    subplot(2,5,4)
    line([min(t), max(t)], [1 1], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t,facilitatory.CAP.AND);
    set(gca, 'XLim', [0 500], 'YLim', [-.5 6], 'FontSize', 12, 'XTick', [])
    ylabel('C_{AND}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    toc
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
end
%% Get ccf
YL = .02;
YH = .1;
XL = .02;
XH = .1;

cdrift = [A YH; A YL; B XH; B XL]; % Note that this line is correct because the first column is the target column
% cdrift = [A YH; A YL; XH B; XL B];
if runCCFf
    for i = 1:numel(pAB)
        for j = 1:size(cdrift,1)
            % Note that when there are two channels that are in conflict, if
            % there is a facilitatory connection, the two conflicting channels
            % will inhibit each other (e.g., -1 * 1 = -1)
            facilitatory.ccdf.ccf(:,j,i) = inh_conf([cdrift(j,:), pAB(i), pBA(i), cA, cB], t);
        end
        facilitatory.ccdf.ccf(facilitatory.ccdf.ccf<0) = 0;
 
        facilitatory.CCF(:,i) =...
            (log(1-facilitatory.ccdf.ccf(:,1,i)) - log(1-facilitatory.ccdf.ccf(:,2,i))) +...
            (log(1-facilitatory.ccdf.ccf(:,3,i)) - log(1-facilitatory.ccdf.ccf(:,4,i)));
    end
    
    %% Plot results
    subplot(2,5,5)
    line([min(t), max(t)], [0 0], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t, facilitatory.CCF);
    set(gca,  'XLim', [0 1500], 'YLim', [-1 5], 'FontSize', 12, 'XTick', [])
    ylabel('CCF(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    legend(h, num2str(pAB'))
    toc
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
end
%%
hA = .05; % High Drift in channel A
hB = .05; % High Drift in channel B
lA = .02; % Low  Drift in channel A
lB = .02; % Low  Drift in channel B

pAB = [0 .25 .5 .75 1];
pBA = pAB;

cA = 10; % Criterion A
cB = 10; % Criterion B

% t = 0:15000; % Time vector to evaluate

drift = [hA, hB; hA lB; lA hB; lA lB];

%% Get cdfs and SICs
if runSICi
    for i = 1:numel(pAB)
        for j = 1:size(drift,1)
            inhibitory.cdf.or(:,j,i) = inh_or([drift(j,:), pAB(i), pBA(i), cA, cB]', t)';
            inhibitory.cdf.and(:,j,i) = inh_and([drift(j,:), pAB(i), pBA(i), cA, cB]', t)';
        end
        inhibitory.SIC.OR(:,i) =...
            (1-inhibitory.cdf.or(:,4,i)) - (1-inhibitory.cdf.or(:,2,i)) -...
            (1-inhibitory.cdf.or(:,3,i)) + (1-inhibitory.cdf.or(:,1,i));
        inhibitory.SIC.AND(:,i) =...
            (1-inhibitory.cdf.and(:,4,i)) - (1-inhibitory.cdf.and(:,2,i)) -...
            (1-inhibitory.cdf.and(:,3,i)) + (1-inhibitory.cdf.and(:,1,i));
    end
    
    %% Plot results
    subplot(2,5,6)
    line([min(t), max(t)], [0 0], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t,inhibitory.SIC.OR);
    set(gca, 'XLim', [0 1500], 'YLim', [-.75 1.25], 'FontSize', 12, 'XTick', [])
    ylabel('SIC_{OR}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
    
    subplot(2,5,7)
    line([min(t), max(t)], [0 0], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t,inhibitory.SIC.AND);
    set(gca, 'XLim', [0 2000], 'YLim', [-.75 1.25], 'FontSize', 12, 'XTick', [])
    ylabel('SIC_{AND}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    toc
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
end

%% Get Capacity
A = .02;
B = .02;
if runCi
    for i = 1:numel(pAB)
        inhibitory.ccdf.OR(:,1,i) = inh_or([A, B, pAB(i), pBA(i), cA, cB]', t)';
        inhibitory.ccdf.OR(:,2,i) = ind_sing([A, cA]', t)';
        inhibitory.ccdf.OR(:,3,i) = ind_sing([B, cB]', t)';
        
        inhibitory.ccdf.AND(:,1,i) = inh_and([A, B, pAB(i), pBA(i), cA, cB]', t)';
        inhibitory.ccdf.AND(:,2,i) = inhibitory.ccdf.OR(:,2,i);
        inhibitory.ccdf.AND(:,3,i) = inhibitory.ccdf.OR(:,3,i);
        
        inhibitory.ccdf.AND(inhibitory.ccdf.AND < 0) = 0;
        
        inhibitory.CAP.OR(:,i) =...
            log(1 - inhibitory.ccdf.OR(:,1,i))./...
            log((1 - inhibitory.ccdf.OR(:,2,i)) .* (1 - inhibitory.ccdf.OR(:,3,i)));      % AB/(A + B);
        
        % There's an issue with the approximation such that the C_AND(t) =
        % 0 prediction when pAB = pBA = 0 doesn't hold. I think this arises
        % from mixing the matrix method with the single channel method
        if i == 1
            inhibitory.CAP.AND(:,i) =...
                log((facilitatory.ccdf.AND(:,2,i)) .* (facilitatory.ccdf.AND(:,3,i)))./...
                log(facilitatory.ccdf.AND(:,1,i));      % (A + B)./AB;
        else
            inhibitory.CAP.AND(:,i) =...
                log((inhibitory.ccdf.AND(:,2,i)) .* (inhibitory.ccdf.AND(:,3,i)))./...
                log(inhibitory.ccdf.AND(:,1,i));      % (A + B)./AB;
        end
    end
    
    %% Plot results
    subplot(2,5,8)
    line([min(t), max(t)], [1 1], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t,inhibitory.CAP.OR);
    set(gca, 'XLim', [0 1500], 'YLim', [-1 6], 'FontSize', 12, 'XTick', [])
    ylabel('C_{OR}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
    
    subplot(2,5,9)
    line([min(t), max(t)], [1 1], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t,inhibitory.CAP.AND);
    set(gca, 'XLim', [0 1200], 'YLim', [-1 6], 'FontSize', 12, 'XTick', [])
    ylabel('C_{AND}(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    toc
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
end

%% Get ccf
YL = .02;
YH = .05;
XL = .02;
XH = .05;

cdrift = [A YH; A YL; B XH; B XL];
% cdrift = [A YH; A YL; XH B; XL B];
if runCCFi
    for i = 1:numel(pAB)
        for j = 1:size(cdrift,1)
            % Note that when there are two channels that are in conflict, if
            % there is an inhibitory connection, the two conflicting channels
            % will facilitate each other (e.g., -(-1) * 1 = 1)
            inhibitory.ccdf.ccf(:,j,i) = fac_conf([cdrift(j,:), pAB(i), pBA(i), cA, cB], t);
        end
        inhibitory.CCF(:,i) =...
            (log(1-inhibitory.ccdf.ccf(:,1,i)) - log(1-inhibitory.ccdf.ccf(:,2,i))) +...
            (log(1-inhibitory.ccdf.ccf(:,3,i)) - log(1-inhibitory.ccdf.ccf(:,4,i)));
    end
    
    %% Plot results
    subplot(2,5,10)
    line([min(t), max(t)], [0 0], 'LineStyle', '--', 'Color', 'k')
    hold on
    h = plot(t,inhibitory.CCF);
    set(gca,  'XLim', [0 1250], 'YLim', [-5 1], 'FontSize', 12, 'XTick', [])
    ylabel('CCF(t)', 'FontSize', 12)
    xlabel('t', 'FontSize', 12)
    box on
    for i = 1:size(colors,1)
        set(h(i), 'Color', colors(i,:), 'LineWidth', 2);
    end
end