%% Figure 10: Exact energy-outage probability (corrected)
% Implements the exact energy-outage expressions of Section 4.1 for Format-I alpha-eta-mu fading.
%
% Important correction:
%   chi_TSR = zeta*rho       (there is NO factor 1/2 for TSR)
%   chi_PSR = zeta*lambda/2
%   chi_IRR = zeta/2
%
% Under zeta = 0.8, rho = 0.3 and lambda = 0.6,
% chi_TSR = chi_PSR = 0.24. Therefore, the TSR and PSR
% energy-outage probabilities must be bit-identical.
%
% The analytical CDF uses the normalized gamma-mixture form of Eq. (18).
% The Monte Carlo generator is valid for arbitrary real mu > 0 in Format I.

clear; close all; clc;
rng(8, 'twister');

%% Common fading and system parameters
alpha = 2;
eta   = 0.5;       % Format-I scattered-power ratio
mu    = 2;
Omega = 1;

zeta   = 0.8;
rho    = 0.3;
lambda = 0.6;
EbarTh = 20;       % E_th/(N0*T)
NsMain = 3;

chiTSR = zeta*rho;
chiPSR = zeta*lambda/2;
chiIRR = zeta/2;
chi    = [chiTSR, chiPSR, chiIRR];

assert(abs(chiTSR-chiPSR) < 10*eps, ...
    'The selected TSR and PSR harvesting coefficients should coincide.');

%% Numerical controls
seriesTol = 1e-12;
nMC       = 1e6;
chunkSize = 1e5;   % limits Monte Carlo memory use

snrDb   = 0:0.1:30;     % smooth analytical curves
snrLin  = 10.^(snrDb/10);
snrDbMC = 0:2:30;       % Monte Carlo marker locations
snrLinMC = 10.^(snrDbMC/10);

%% Exact analytical energy-outage probabilities
xTSR = EbarTh./(chiTSR*snrLin);
xPSR = EbarTh./(chiPSR*snrLin);
xIRR = EbarTh./(chiIRR*snrLin);

[Ftsr, Kused, seriesResidual] = alphaEtaMuCdfFormatI( ...
    xTSR, alpha, eta, mu, Omega, seriesTol);
Fpsr = alphaEtaMuCdfFormatI( ...
    xPSR, alpha, eta, mu, Omega, seriesTol);
Firr = alphaEtaMuCdfFormatI( ...
    xIRR, alpha, eta, mu, Omega, seriesTol);

PE_TSR = Ftsr.^NsMain;
PE_PSR = Fpsr.^NsMain;
PE_IRR = Firr.^NsMain;

assert(max(abs(PE_TSR-PE_PSR)) < 1e-14, ...
    'TSR and PSR analytical energy-outage curves must coincide.');

%% Inset: source-side candidate-position selection
NsInset = [1, 3, 5];
Finset = alphaEtaMuCdfFormatI( ...
    EbarTh./(chiPSR*snrLin), alpha, eta, mu, Omega, seriesTol);
PE_inset = Finset(:).^(NsInset);  % each column corresponds to one Ns

%% Monte Carlo validation using common random numbers
% Common channel realizations make the equality of TSR and PSR visible
% without artificial differences caused by independent simulation runs.
[PE_MC_main, PE_MC_inset] = monteCarloEnergyOutage( ...
    nMC, chunkSize, max(NsInset), NsMain, NsInset, ...
    alpha, eta, mu, Omega, EbarTh, chi, snrLinMC);

assert(isequal(PE_MC_main(1,:), PE_MC_main(2,:)), ...
    'TSR and PSR Monte Carlo estimates should be exactly identical.');

% Analytical values at the Monte Carlo marker locations (diagnostics only)
PE_exact_MC = zeros(3, numel(snrLinMC));
for p = 1:3
    Fp = alphaEtaMuCdfFormatI(EbarTh./(chi(p)*snrLinMC), ...
        alpha, eta, mu, Omega, seriesTol);
    PE_exact_MC(p,:) = Fp.^NsMain;
end

reliableMask = nMC*PE_exact_MC >= 20; % avoid judging extremely rare-event MC
absError = abs(PE_MC_main-PE_exact_MC);
fprintf('Figure 8 diagnostics\n');
fprintf('  chi_TSR = %.4f, chi_PSR = %.4f, chi_IRR = %.4f\n', ...
    chiTSR, chiPSR, chiIRR);
fprintf('  CDF truncation order K = %d\n', Kused);
fprintf('  Remaining mixture-weight residual = %.3e\n', seriesResidual);
fprintf('  Maximum MC absolute error where >=20 events are expected = %.3e\n', ...
    max(absError(reliableMask)));

%% Plot
% The following settings match the publication style used in Figures 2-7.
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 3.5 3.10], ...
    'PaperPositionMode','auto');

ax = axes(fig);
hold(ax,'on');
box(ax,'on');
grid(ax,'on');

% Protocol colours are kept identical to the earlier comparison figure:
% IRR = blue, TSR = orange and PSR = green.
cIRR = [0.0000 0.4470 0.7410];
cTSR = [0.8500 0.3250 0.0980];
cPSR = [0.4660 0.6740 0.1880];

% Plot PSR first with a slightly wider line and overlay TSR.  This makes
% their exact coincidence visible without shifting either curve.
hPSR = semilogy(ax,snrDb,PE_PSR,'-.', ...
    'Color',cPSR, ...
    'LineWidth',2.0, ...
    'DisplayName','PSR, analytical');

hTSR = semilogy(ax,snrDb,PE_TSR,'--', ...
    'Color',cTSR, ...
    'LineWidth',1.2, ...
    'DisplayName','TSR, analytical');

hIRR = semilogy(ax,snrDb,PE_IRR,'-', ...
    'Color',cIRR, ...
    'LineWidth',1.2, ...
    'DisplayName','IRR, analytical');

% Omit zero-event estimates on the logarithmic scale.  All simulation
% points use the open black circles employed in Figures 2-6.
mcCommon = zeroToNaN(PE_MC_main(1,:));
mcIRR    = zeroToNaN(PE_MC_main(3,:));

hSimulation = semilogy(ax,snrDbMC,mcCommon,'ko', ...
    'LineStyle','none', ...
    'MarkerFaceColor','none', ...
    'MarkerSize',5, ...
    'LineWidth',0.8, ...
    'DisplayName','Simulation');

semilogy(ax,snrDbMC,mcIRR,'ko', ...
    'LineStyle','none', ...
    'MarkerFaceColor','none', ...
    'MarkerSize',5, ...
    'LineWidth',0.8, ...
    'HandleVisibility','off');

xlabel(ax,'Average source transmit SNR, $\bar{\gamma}$ (dB)', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

ylabel(ax,'Energy-outage probability, $P_{\mathrm{E},\nu}$', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

xlim(ax,[0 30]);
xticks(ax,0:5:30);
ylim(ax,[1e-6 1]);
yticks(ax,10.^(-6:1:0));

set(ax, ...
    'YScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',0.8, ...
    'TickDir','in', ...
    'XMinorGrid','off', ...
    'YMinorGrid','on', ...
    'GridAlpha',0.22, ...
    'MinorGridAlpha',0.12, ...
    'Layer','top');

lgd = legend(ax,[hIRR,hTSR,hPSR,hSimulation], ...
    {'IRR, analytical','TSR, analytical', ...
     'PSR, analytical','Simulation'}, ...
    'Location','southwest', ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'Interpreter','latex', ...
    'Box','on');

lgd.LineWidth = 1;
lgd.ItemTokenSize = [11 11];

%% Inset axes: Ns = 1, 3 and 5 for PSR
axIn = axes(fig,'Position',[0.28 0.24 0.41 0.30]);
hold(axIn,'on');
box(axIn,'on');
grid(axIn,'on');

insetColors = [cIRR; cTSR; cPSR];
insetStyles = {'-','--','-.'};
hInset = gobjects(1,numel(NsInset));
for j = 1:numel(NsInset)
    hInset(j) = semilogy(axIn,snrDb,PE_inset(:,j), ...
        insetStyles{j}, ...
        'Color',insetColors(j,:), ...
        'LineWidth',1.0);

    semilogy(axIn,snrDbMC,zeroToNaN(PE_MC_inset(j,:)),'ko', ...
        'LineStyle','none', ...
        'MarkerFaceColor','none', ...
        'MarkerSize',3, ...
        'LineWidth',0.6, ...
        'HandleVisibility','off');
end

xlim(axIn,[12 24]);
xticks(axIn,12:4:24);
ylim(axIn,[1e-7 1]);
yticks(axIn,10.^(-6:2:0));

xlabel(axIn,'SNR (dB)', ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'Interpreter','latex');

ylabel(axIn,'$P_{\mathrm{E},\mathrm{PSR}}$', ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'Interpreter','latex');

set(axIn, ...
    'YScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'LineWidth',0.65, ...
    'TickDir','in', ...
    'XMinorGrid','off', ...
    'YMinorGrid','on', ...
    'GridAlpha',0.22, ...
    'MinorGridAlpha',0.12, ...
    'Layer','top');

lgdInset = legend(axIn,hInset, ...
    {'$N_s=1$','$N_s=3$','$N_s=5$'}, ...
    'Interpreter','latex', ...
    'Location','southwest', ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'Box','off');

lgdInset.ItemTokenSize = [9 9];

%% Export publication files
outputDir = fullfile(pwd, 'Fig10_output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if exist('exportgraphics','file') == 2
    exportgraphics(fig, fullfile(outputDir, 'Fig10_Energy_Outage.pdf'), ...
        'ContentType','vector');
    exportgraphics(fig, fullfile(outputDir, 'Fig10_Energy_Outage.png'), ...
        'Resolution',600);
else
    % Fallback for MATLAB releases older than R2020a.
    print(fig, fullfile(outputDir, 'Fig10_Energy_Outage.pdf'), ...
        '-dpdf','-painters');
    print(fig, fullfile(outputDir, 'Fig10_Energy_Outage.png'), ...
        '-dpng','-r600');
end
savefig(fig, fullfile(outputDir, 'Fig10_Energy_Outage.fig'));

%% Local functions
function [F, K, residual] = alphaEtaMuCdfFormatI( ...
    z, alpha, eta, mu, Omega, tol)
% Parent alpha-eta-mu power-gain CDF, evaluated as a positive mixture of
% regularized lower incomplete gamma functions.
%
% MATLAB gammainc(x,a,'lower') returns gamma(a,x)/Gamma(a), so the
% unregularized Gamma(a) factor is absorbed into the mixture weights.

    validateattributes(z,     {'numeric'}, {'real','nonnegative'});
    validateattributes(alpha, {'numeric'}, {'scalar','real','positive'});
    validateattributes(eta,   {'numeric'}, {'scalar','real','positive'});
    validateattributes(mu,    {'numeric'}, {'scalar','real','positive'});
    validateattributes(Omega, {'numeric'}, {'scalar','real','positive'});

    h = (1+eta)^2/(4*eta);
    H = (1-eta^2)/(4*eta);

    % Normalized gamma-mixture weights:
    % w_0 = h^(-mu),
    % w_(k+1)/w_k = (H/h)^2 * (mu+k)/(k+1).
    maxK = 100000;
    weights = zeros(1,64);
    weights(1) = h^(-mu);
    sumWeights = weights(1);
    K = 0;

    while (1-sumWeights) > tol
        if K >= maxK
            error('alphaEtaMuCdfFormatI:NoConvergence', ...
                'Series did not meet the requested tolerance.');
        end
        if K+2 > numel(weights)
            weights = [weights, zeros(1,numel(weights))]; %#ok<AGROW>
        end
        ratio = (H/h)^2 * (mu+K)/(K+1);
        weights(K+2) = weights(K+1)*ratio;
        K = K+1;
        sumWeights = sumWeights + weights(K+1);
    end

    weights = weights(1:K+1);
    residual = max(0,1-sumWeights);

    x = 2*mu*h*(z./Omega).^(alpha/2);
    F = zeros(size(z));
    for k = 0:K
        shape = 2*(mu+k);
        F = F + weights(k+1)*gammainc(x,shape,'lower');
    end
    F = min(max(F,0),1);
end

function Z = sampleAlphaEtaMuFormatI(nRows,nCols,alpha,eta,mu,Omega)
% Exact Format-I sampler for arbitrary real mu > 0.
% If S = Gx+Gy, then Z = S^(2/alpha), where the independent gamma
% components have shape mu and the scales shown below.

    sigmaY2 = Omega^(alpha/2)/(2*mu*(1+eta));
    sigmaX2 = eta*sigmaY2;

    Gx = (2*sigmaX2)*randg(mu,nRows,nCols);
    Gy = (2*sigmaY2)*randg(mu,nRows,nCols);
    Z = (Gx+Gy).^(2/alpha);
end

function [PEmain,PEinset] = monteCarloEnergyOutage( ...
    nMC,chunkSize,NsMax,NsMain,NsInset, ...
    alpha,eta,mu,Omega,EbarTh,chi,snrLin)
% Chunked Monte Carlo evaluation using common random numbers.

    nProtocol = numel(chi);
    nSNR = numel(snrLin);
    countMain  = zeros(nProtocol,nSNR);
    countInset = zeros(numel(NsInset),nSNR);

    thresholds = EbarTh./(chi(:)*snrLin(:).');
    insetThresholds = EbarTh./(chi(2)*snrLin(:).');

    nCompleted = 0;
    while nCompleted < nMC
        nNow = min(chunkSize,nMC-nCompleted);
        Z = sampleAlphaEtaMuFormatI( ...
            nNow,NsMax,alpha,eta,mu,Omega);
        Zmax = cummax(Z,2);

        Xmain = Zmax(:,NsMain);
        for p = 1:nProtocol
            countMain(p,:) = countMain(p,:) + ...
                sum(bsxfun(@lt,Xmain,thresholds(p,:)),1);
        end

        for j = 1:numel(NsInset)
            Xinset = Zmax(:,NsInset(j));
            countInset(j,:) = countInset(j,:) + ...
                sum(bsxfun(@lt,Xinset,insetThresholds),1);
        end
        nCompleted = nCompleted+nNow;
    end

    PEmain  = countMain/nMC;
    PEinset = countInset/nMC;
end

function y = zeroToNaN(y)
% Zero-event estimates cannot be displayed on a logarithmic axis.
    y(y==0) = NaN;
end
