%% Figure 12: Optimized throughput versus normalized energy demand
% This script reproduces the optimized-throughput versus energy-demand study in the revised manuscript.
%
% Main graph:
%   - exact optimized throughput for IRR, PSR, and TSR;
%   - Monte Carlo markers, following the style of Figs. 1-9;
%   - energy-constraint activation and feasibility limits.
%
% The output is sized for one column of a two-column paper.

clear; close all; clc;
rng(10,'twister');

%% System and channel parameters
% Homogeneous Format-I alpha-eta-mu fading on the two hops.
alpha1 = 2;   eta1 = 0.5;   mu1 = 2;   Omega1 = 1;
alpha2 = 2;   eta2 = 0.5;   mu2 = 2;   Omega2 = 1;

Ns = 3;
Nd = 3;
R = 2;                              % bit/s/Hz over the complete block
zeta = 0.8;
gammaBarDb = 20;
gammaBar = 10^(gammaBarDb/10);      % P_S/N_0
epsilonE = 1e-2;

seriesTol = 1e-13;
dist1 = makeAlphaEtaMuFormatI(alpha1,eta1,mu1,Omega1,seriesTol);
dist2 = makeAlphaEtaMuFormatI(alpha2,eta2,mu2,Omega2,seriesTol);

%% Energy-feasibility constants
qNs = parentQuantile(epsilonE^(1/Ns),dist1);
energyScale = zeta*gammaBar*qNs;

% TSR: rho_min = Ebar_th/energyScale.
% PSR: lambda_min = 2*Ebar_th/energyScale.
% IRR has the same harvesting coefficient as PSR at lambda = 1.
EmaxTSR = energyScale;
EmaxPSR = energyScale/2;
EmaxIRR = EmaxPSR;

%% Unconstrained TSR and PSR optima
fprintf('\nFigure 12: optimized throughput versus energy demand\n');
fprintf('  q_Ns = %.10f\n',qNs);

rhoCoarse = linspace(0.01,0.60,91);
lambdaCoarse = linspace(0.01,0.99,111);
TtsrCoarse = zeros(size(rhoCoarse));
TpsrCoarse = zeros(size(lambdaCoarse));

fprintf('  Finding the unconstrained TSR optimum ...\n');
for k = 1:numel(rhoCoarse)
    TtsrCoarse(k) = exactThroughput('TSR',rhoCoarse(k),R,zeta, ...
        gammaBar,Ns,Nd,dist1,dist2);
end
[rhoStar,TtsrStar] = refineMaximum('TSR',rhoCoarse,TtsrCoarse, ...
    R,zeta,gammaBar,Ns,Nd,dist1,dist2);

fprintf('  Finding the unconstrained PSR optimum ...\n');
for k = 1:numel(lambdaCoarse)
    TpsrCoarse(k) = exactThroughput('PSR',lambdaCoarse(k),R,zeta, ...
        gammaBar,Ns,Nd,dist1,dist2);
end
[lambdaStar,TpsrStar] = refineMaximum('PSR',lambdaCoarse, ...
    TpsrCoarse,R,zeta,gammaBar,Ns,Nd,dist1,dist2);

% IRR has no adjustable protocol factor.
TirrExact = exactThroughput('IRR',NaN,R,zeta,gammaBar, ...
    Ns,Nd,dist1,dist2);

% Energy demands at which the constraints first become active.
EbindTSR = rhoStar*energyScale;
EbindPSR = lambdaStar*energyScale/2;

fprintf('\n  Unconstrained optima and transition points\n');
fprintf('    TSR: rho* = %.7f, T* = %.7f, constraint active at Ebar_th = %.5f\n', ...
    rhoStar,TtsrStar,EbindTSR);
fprintf('    PSR: lambda* = %.7f, T* = %.7f, constraint active at Ebar_th = %.5f\n', ...
    lambdaStar,TpsrStar,EbindPSR);
fprintf('    IRR: T = %.7f\n',TirrExact);
fprintf('    PSR/IRR feasibility limit = %.5f\n',EmaxPSR);
fprintf('    TSR feasibility limit = %.5f\n',EmaxTSR);

%% Exact optimized-throughput curves
% The plotted range stops at 25 because the TSR throughput is already
% negligible there. Its formal feasibility limit is printed above.
EplotMax = 25;
EnearPSRLimit = EmaxPSR*(1-1e-7);
EbarGrid = unique([linspace(0,EplotMax,151),EbindTSR,EbindPSR, ...
                   EnearPSRLimit,EmaxPSR]);

rhoOpt = nan(size(EbarGrid));
lambdaOpt = nan(size(EbarGrid));
TtsrOpt = nan(size(EbarGrid));
TpsrOpt = nan(size(EbarGrid));
TirrOpt = nan(size(EbarGrid));

fprintf('\n  Evaluating optimized energy-demand curves ...\n');
for k = 1:numel(EbarGrid)
    Ebar = EbarGrid(k);

    rhoMin = Ebar/energyScale;
    if rhoMin < 1
        if rhoMin <= rhoStar
            rhoOpt(k) = rhoStar;
            TtsrOpt(k) = TtsrStar;
        else
            rhoOpt(k) = rhoMin;
            TtsrOpt(k) = exactThroughput('TSR',rhoOpt(k),R,zeta, ...
                gammaBar,Ns,Nd,dist1,dist2);
        end
    end

    lambdaMin = 2*Ebar/energyScale;
    if lambdaMin < 1
        if lambdaMin <= lambdaStar
            lambdaOpt(k) = lambdaStar;
            TpsrOpt(k) = TpsrStar;
        else
            lambdaOpt(k) = lambdaMin;
            TpsrOpt(k) = exactThroughput('PSR',lambdaOpt(k),R,zeta, ...
                gammaBar,Ns,Nd,dist1,dist2);
        end
    end

    if Ebar <= EmaxIRR
        TirrOpt(k) = TirrExact;
    end
end

%% Monte Carlo validation
% The optimized protocol factors are taken from the analytical energy
% constraint, while the information-outage probability is simulated.
EbarMC = [0,4,8,10,12,14,16,18,20,22,23];
rhoMC = max(rhoStar,EbarMC/energyScale);
lambdaMC = max(lambdaStar,2*EbarMC/energyScale);
rhoMC(rhoMC >= 1) = NaN;
lambdaMC(lambdaMC >= 1) = NaN;

nMC = 1e6;
chunkSize = 1e5;
fprintf('\n  Running Monte Carlo validation with %d realizations ...\n',nMC);
[TtsrMC,TpsrMC,TirrMC] = monteCarloOptimizedThroughput( ...
    rhoMC,lambdaMC,EbarMC,EmaxIRR,nMC,chunkSize,R,zeta,gammaBar, ...
    Ns,Nd,alpha1,eta1,mu1,Omega1,alpha2,eta2,mu2,Omega2);

% Exact values at the Monte Carlo marker positions.
TtsrAtMC = nan(size(EbarMC));
TpsrAtMC = nan(size(EbarMC));
TirrAtMC = nan(size(EbarMC));
for k = 1:numel(EbarMC)
    if isfinite(rhoMC(k))
        TtsrAtMC(k) = exactThroughput('TSR',rhoMC(k),R,zeta, ...
            gammaBar,Ns,Nd,dist1,dist2);
    end
    if isfinite(lambdaMC(k))
        TpsrAtMC(k) = exactThroughput('PSR',lambdaMC(k),R,zeta, ...
            gammaBar,Ns,Nd,dist1,dist2);
    end
    if EbarMC(k) <= EmaxIRR
        TirrAtMC(k) = TirrExact;
    end
end

fprintf('    Maximum |MC-exact|, TSR = %.3e\n', ...
    maxFinite(abs(TtsrMC-TtsrAtMC)));
fprintf('    Maximum |MC-exact|, PSR = %.3e\n', ...
    maxFinite(abs(TpsrMC-TpsrAtMC)));
fprintf('    Maximum |MC-exact|, IRR = %.3e\n', ...
    maxFinite(abs(TirrMC-TirrAtMC)));

%% Plot: one-column journal layout
% The compact legend is placed in the unused lower-left region so that
% the binding transitions and feasibility limits remain visible.
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

% Protocol colours and line styles match Figures 2 and 8:
% IRR: blue solid; TSR: orange dashed; PSR: green dash-dot.
blue = [0.0000,0.4470,0.7410];
orange = [0.8500,0.3250,0.0980];
green = [0.4660,0.6740,0.1880];

hIRR = plot(ax,EbarGrid,TirrOpt,'-', ...
    'Color',blue,'LineWidth',1.2);
hTSR = plot(ax,EbarGrid,TtsrOpt,'--', ...
    'Color',orange,'LineWidth',1.2);
hPSR = plot(ax,EbarGrid,TpsrOpt,'-.', ...
    'Color',green,'LineWidth',1.2);

% Open black circles match the simulation style of Figures 2-9.
hMC = plot(ax,EbarMC,TirrMC,'ko','LineStyle','none', ...
    'MarkerFaceColor','none','MarkerSize',5,'LineWidth',0.8);
plot(ax,EbarMC,TtsrMC,'ko','LineStyle','none', ...
    'MarkerFaceColor','none','MarkerSize',5,'LineWidth',0.8, ...
    'HandleVisibility','off');
plot(ax,EbarMC,TpsrMC,'ko','LineStyle','none', ...
    'MarkerFaceColor','none','MarkerSize',5,'LineWidth',0.8, ...
    'HandleVisibility','off');

% As in Figure 9, identical filled black stars mark the key transition
% points.  Here they show where each energy constraint becomes active.
hBinding = plot(ax,EbindTSR,TtsrStar,'p', ...
    'Color','k', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',9, ...
    'LineWidth',0.8, ...
    'LineStyle','none');

plot(ax,EbindPSR,TpsrStar,'p', ...
    'Color','k', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',9, ...
    'LineWidth',0.8, ...
    'LineStyle','none', ...
    'HandleVisibility','off');

% PSR becomes inadmissible and IRR reaches its feasibility limit at the
% same normalized energy demand. IRR includes the endpoint; PSR does not.
plot(ax,[EmaxPSR,EmaxPSR],[0,2.04],':','Color',[0.35,0.35,0.35], ...
    'LineWidth',0.8,'HandleVisibility','off');
plot(ax,EmaxPSR,0,'d','Color',green,'MarkerFaceColor','w', ...
    'MarkerSize',5,'LineWidth',0.8,'HandleVisibility','off');
text(ax,EmaxPSR-0.25,0.08,'PSR/IRR limit', ...
    'Rotation',90,'HorizontalAlignment','left', ...
    'Color',[0.30,0.30,0.30], ...
    'FontName','Times New Roman','FontSize',8);

xlabel(ax,'Normalized energy threshold, $\bar{E}_{\mathrm{th}}$', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');
ylabel(ax,'Optimized throughput, $\mathcal{T}_{\nu}^{\star}$ (bit/s/Hz)', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');
xlim(ax,[0,EplotMax]);
xticks(ax,0:5:EplotMax);
ylim(ax,[0,2.08]);

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',0.8, ...
    'TickDir','in', ...
    'XMinorGrid','off', ...
    'YMinorGrid','on', ...
    'GridAlpha',0.22, ...
    'MinorGridAlpha',0.12, ...
    'Layer','top');

% Compact legend inside the unused lower-left part of the axes.
leg = legend(ax,[hIRR,hTSR,hPSR,hMC,hBinding], ...
    {'IRR','TSR','PSR','Simulation','Binding point'}, ...
    'Location','southwest', ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'Box','on', ...
    'Color','w', ...
    'EdgeColor','k', ...
    'Interpreter','tex');

leg.LineWidth = 1;
leg.ItemTokenSize = [11 11];

%% Export publication files and numerical tables
outputDir = fullfile(pwd,'Fig12_output');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

if exist('exportgraphics','file') == 2
    exportgraphics(fig,fullfile(outputDir, ...
        'Fig12_optimized_energy_demand.pdf'),'ContentType','vector');
    exportgraphics(fig,fullfile(outputDir, ...
        'Fig12_optimized_energy_demand.png'),'Resolution',600);
else
    print(fig,fullfile(outputDir, ...
        'Fig12_optimized_energy_demand.pdf'),'-dpdf','-painters');
    print(fig,fullfile(outputDir, ...
        'Fig12_optimized_energy_demand.png'),'-dpng','-r600');
end
savefig(fig,fullfile(outputDir,'Fig12_optimized_energy_demand.fig'));

curveTable = table(EbarGrid(:),rhoOpt(:),TtsrOpt(:), ...
    lambdaOpt(:),TpsrOpt(:),TirrOpt(:), ...
    'VariableNames',{'Ebar_th','rho_opt','T_TSR_opt', ...
                     'lambda_opt','T_PSR_opt','T_IRR'});
writetable(curveTable,fullfile(outputDir,'Fig12_curves.csv'));

mcTable = table(EbarMC(:),TtsrMC(:),TpsrMC(:),TirrMC(:), ...
    'VariableNames',{'Ebar_th','T_TSR_MC','T_PSR_MC','T_IRR_MC'});
writetable(mcTable,fullfile(outputDir,'Fig12_monte_carlo.csv'));

keyTable = table(qNs,rhoStar,lambdaStar,EbindTSR,EbindPSR, ...
    EmaxTSR,EmaxPSR,TtsrStar,TpsrStar,TirrExact);
writetable(keyTable,fullfile(outputDir,'Fig12_key_values.csv'));

%% Local functions
function value = maxFinite(x)
% Maximum after excluding NaN and Inf values (older MATLAB compatible).
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = max(x);
    end
end

function dist = makeAlphaEtaMuFormatI(alpha,eta,mu,Omega,tol)
% Positive gamma-mixture representation of the Format-I parent CDF.
    h = (1+eta)^2/(4*eta);
    H = (1-eta^2)/(4*eta);

    weights = zeros(1,64);
    weights(1) = h^(-mu);
    sumWeights = weights(1);
    k = 0;
    maxK = 100000;

    while (1-sumWeights) > tol
        if k >= maxK
            error('makeAlphaEtaMuFormatI:NoConvergence', ...
                'The alpha-eta-mu CDF series did not converge.');
        end
        if k+2 > numel(weights)
            weights = [weights,zeros(1,numel(weights))]; %#ok<AGROW>
        end
        ratio = (H/h)^2*(mu+k)/(k+1);
        weights(k+2) = weights(k+1)*ratio;
        k = k+1;
        sumWeights = sumWeights+weights(k+1);
    end

    weights = weights(1:k+1);
    residual = max(0,1-sumWeights);
    weights = weights/sum(weights);

    dist.alpha = alpha;
    dist.eta = eta;
    dist.mu = mu;
    dist.Omega = Omega;
    dist.h = h;
    dist.H = H;
    dist.weights = weights;
    dist.shapes = 2*(mu+(0:k));
    dist.rate = 2*mu*h/Omega^(alpha/2);
    dist.K = k;
    dist.seriesResidual = residual;
end

function F = parentCdf(z,dist)
% Parent channel-power CDF in Eq. (18).
    F = zeros(size(z));
    positive = z > 0;
    if any(positive(:))
        x = dist.rate*z(positive).^(dist.alpha/2);
        temp = zeros(size(x));
        for k = 1:numel(dist.weights)
            temp = temp+dist.weights(k)* ...
                gammainc(x,dist.shapes(k),'lower');
        end
        F(positive) = temp;
    end
    F(isinf(z) & z>0) = 1;
    F = min(max(F,0),1);
end

function f = parentPdf(z,dist)
% Parent channel-power PDF from the same gamma mixture.
    f = zeros(size(z));
    positive = isfinite(z) & z > 0;
    if ~any(positive(:))
        return;
    end

    zp = z(positive);
    x = dist.rate*zp.^(dist.alpha/2);
    dxDz = dist.rate*(dist.alpha/2)*zp.^(dist.alpha/2-1);
    temp = zeros(size(x));
    logx = log(x);
    for k = 1:numel(dist.weights)
        shape = dist.shapes(k);
        gammaPdf = exp((shape-1).*logx-x-gammaln(shape));
        temp = temp+dist.weights(k)*gammaPdf.*dxDz;
    end
    f(positive) = temp;
    f = max(f,0);
end

function F = selectedCdf(z,N,dist)
% CDF of the maximum of N independent candidate-position gains.
    F = parentCdf(z,dist).^N;
end

function f = selectedPdf(z,N,dist)
% PDF of the maximum of N independent candidate-position gains.
    Fz = parentCdf(z,dist);
    f = N*parentPdf(z,dist).*Fz.^(N-1);
end

function q = parentQuantile(p,dist)
% Generalized inverse of the parent CDF, computed by bisection.
    if p <= 0
        q = 0;
        return;
    elseif p >= 1
        q = Inf;
        return;
    end

    lower = 0;
    upper = max(1,dist.Omega);
    while parentCdf(upper,dist) < p
        upper = 2*upper;
        if upper > 1e12*max(1,dist.Omega)
            error('parentQuantile:NoBracket','Unable to bracket quantile.');
        end
    end

    for iter = 1:90
        middle = (lower+upper)/2;
        if parentCdf(middle,dist) < p
            lower = middle;
        else
            upper = middle;
        end
    end
    q = (lower+upper)/2;
end

function [a,b,gammaTh,valid] = protocolParameters( ...
    protocol,theta,R,zeta,gammaBar)
% Protocol substitutions for TSR, PSR, and IRR.
    a = NaN; b = NaN; gammaTh = Inf; valid = true;

    switch upper(protocol)
        case 'TSR'
            valid = isfinite(theta) && theta > 0 && theta < 1;
            if ~valid
                return;
            end
            exponent = 2*R/(1-theta);
            if exponent >= log2(realmax)
                valid = false;
                return;
            end
            a = gammaBar;
            kappa = 2*zeta*theta/(1-theta);
            b = kappa*gammaBar;
            gammaTh = 2^exponent-1;

        case 'PSR'
            valid = isfinite(theta) && theta > 0 && theta < 1;
            if ~valid
                return;
            end
            a = (1-theta)*gammaBar;
            b = zeta*theta*gammaBar;
            gammaTh = 2^(2*R)-1;

        case 'IRR'
            a = gammaBar;
            b = zeta*gammaBar;
            gammaTh = 2^(2*R)-1;

        otherwise
            error('protocolParameters:UnknownProtocol', ...
                'Protocol must be TSR, PSR, or IRR.');
    end
end

function T = exactThroughput(protocol,theta,R,zeta,gammaBar, ...
    Ns,Nd,dist1,dist2)
% Exact delay-limited throughput from the one-dimensional outage representation.
    [a,b,gammaTh,valid] = protocolParameters( ...
        protocol,theta,R,zeta,gammaBar);
    if ~valid || ~isfinite(gammaTh)
        T = 0;
        return;
    end

    xTh = gammaTh/a;
    L1 = selectedCdf(xTh,Ns,dist1);
    if L1 >= 1-1e-13
        T = 0;
        return;
    end

    integrand = @(t) transformedOutageIntegrand( ...
        t,xTh,a,b,gammaTh,Ns,Nd,dist1,dist2);
    L2 = integral(integrand,0,1,'RelTol',2e-8,'AbsTol',1e-11);
    Pout = min(max(L1+L2,0),1);
    T = R*(1-Pout);
end

function value = transformedOutageIntegrand( ...
    t,xTh,a,b,gammaTh,Ns,Nd,dist1,dist2)
% Map x in [x_th,infinity) to t in [0,1) by x=x_th+t/(1-t).
    value = zeros(size(t));
    mask = t >= 0 & t < 1;
    if ~any(mask(:))
        return;
    end

    tm = t(mask);
    x = xTh+tm./(1-tm);
    denominator = b*x.*(a*x-gammaTh);
    G = ones(size(x));
    regular = denominator > 0;
    if any(regular(:))
        yRequired = gammaTh*(a*x(regular)+1)./denominator(regular);
        G(regular) = selectedCdf(yRequired,Nd,dist2);
    end

    fX = selectedPdf(x,Ns,dist1);
    value(mask) = G.*fX./(1-tm).^2;
    value(~isfinite(value)) = 0;
end

function [thetaStar,TStar] = refineMaximum(protocol,grid,Tgrid, ...
    R,zeta,gammaBar,Ns,Nd,dist1,dist2)
% Refine the maximum within neighboring coarse-grid points.
    [~,index] = max(Tgrid);
    indexLower = max(1,index-1);
    indexUpper = min(numel(grid),index+1);
    lower = max(grid(indexLower),1e-6);
    upper = min(grid(indexUpper),1-1e-6);
    objective = @(theta) -exactThroughput(protocol,theta,R,zeta, ...
        gammaBar,Ns,Nd,dist1,dist2);
    options = optimset('TolX',1e-8,'Display','off');
    thetaStar = fminbnd(objective,lower,upper,options);
    TStar = -objective(thetaStar);
end

function [TtsrMC,TpsrMC,TirrMC] = monteCarloOptimizedThroughput( ...
    rhoValues,lambdaValues,EbarValues,EmaxIRR,nMC,chunkSize, ...
    R,zeta,gammaBar,Ns,Nd,alpha1,eta1,mu1,Omega1, ...
    alpha2,eta2,mu2,Omega2)
% Chunked Monte Carlo simulation of the end-to-end AF SNR.
    successTSR = zeros(size(rhoValues));
    successPSR = zeros(size(lambdaValues));
    successIRR = 0;
    completed = 0;

    gammaThFixed = 2^(2*R)-1;
    while completed < nMC
        nNow = min(chunkSize,nMC-completed);
        Z1 = sampleAlphaEtaMuFormatI( ...
            nNow,Ns,alpha1,eta1,mu1,Omega1);
        Z2 = sampleAlphaEtaMuFormatI( ...
            nNow,Nd,alpha2,eta2,mu2,Omega2);
        X = max(Z1,[],2);
        Y = max(Z2,[],2);

        for k = 1:numel(rhoValues)
            rho = rhoValues(k);
            if isfinite(rho)
                a = gammaBar;
                b = (2*zeta*rho/(1-rho))*gammaBar;
                gammaTh = 2^(2*R/(1-rho))-1;
                gammaAF = (a*b*X.^2.*Y)./(a*X+b*X.*Y+1);
                successTSR(k) = successTSR(k)+sum(gammaAF >= gammaTh);
            end
        end

        for k = 1:numel(lambdaValues)
            lambda = lambdaValues(k);
            if isfinite(lambda)
                a = (1-lambda)*gammaBar;
                b = zeta*lambda*gammaBar;
                gammaAF = (a*b*X.^2.*Y)./(a*X+b*X.*Y+1);
                successPSR(k) = successPSR(k)+ ...
                    sum(gammaAF >= gammaThFixed);
            end
        end

        a = gammaBar;
        b = zeta*gammaBar;
        gammaAF = (a*b*X.^2.*Y)./(a*X+b*X.*Y+1);
        successIRR = successIRR+sum(gammaAF >= gammaThFixed);
        completed = completed+nNow;
    end

    TtsrMC = R*successTSR/nMC;
    TpsrMC = R*successPSR/nMC;
    TirrValue = R*successIRR/nMC;
    TirrMC = TirrValue*ones(size(EbarValues));

    TtsrMC(~isfinite(rhoValues)) = NaN;
    TpsrMC(~isfinite(lambdaValues)) = NaN;
    TirrMC(EbarValues > EmaxIRR) = NaN;
end

function Z = sampleAlphaEtaMuFormatI(nRows,nCols,alpha,eta,mu,Omega)
% Exact Format-I sampler for arbitrary real mu > 0.
    sigmaY2 = Omega^(alpha/2)/(2*mu*(1+eta));
    sigmaX2 = eta*sigmaY2;
    Gx = (2*sigmaX2)*randg(mu,nRows,nCols);
    Gy = (2*sigmaY2)*randg(mu,nRows,nCols);
    Z = (Gx+Gy).^(2/alpha);
end
