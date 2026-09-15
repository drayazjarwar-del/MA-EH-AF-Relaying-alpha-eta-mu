%% Figure 13: Joint optimization versus fixed-design baselines
% This script reproduces the joint PSR design problem in the revised manuscript.
%
% Four strategies are compared for a total position budget NT = 8:
%   1) Joint: optimize both Ns (and Nd=NT-Ns) and lambda.
%   2) Allocation-only: optimize Ns with lambda fixed at 0.6.
%   3) Protocol-only: optimize lambda with Ns=Nd=4.
%   4) Fully fixed: use Ns=Nd=4 and lambda=0.6.
%
% Analytical lines use the exact one-dimensional outage representation. Open black
% circles are Monte Carlo simulations, following the style of Figs. 1-10.
% Infeasible designs are terminated rather than plotted as zero.

clear; close all; clc;
rng(11,'twister');

%% Runtime environment (for reproducibility/runtime reporting)
% These lines do not change the numerical calculations. They only record
% the software/hardware environment used for the representative runtime.
matlabVersionText = version;
matlabReleaseText = version('-release');
architectureText = computer('arch');
cpuText = getenv('PROCESSOR_IDENTIFIER');
if isempty(cpuText)
    cpuText = 'Not reported by operating system';
end

fprintf('\nRuntime environment\n');
fprintf('  MATLAB version: %s\n',matlabVersionText);
fprintf('  MATLAB release: %s\n',matlabReleaseText);
fprintf('  Architecture:   %s\n',architectureText);
fprintf('  CPU:            %s\n',cpuText);

%% System and channel parameters
% Homogeneous Format-I alpha-eta-mu fading on the two hops.
alpha1 = 2;   eta1 = 0.5;   mu1 = 2;   Omega1 = 1;
alpha2 = 2;   eta2 = 0.5;   mu2 = 2;   Omega2 = 1;

NT = 8;
R = 2;                              % bit/s/Hz over the complete block
zeta = 0.8;
gammaBarDb = 20;
gammaBar = 10^(gammaBarDb/10);
EbarTh = 12;                        % E_th/(N_0 T)
lambdaFixed = 0.6;

seriesTol = 1e-13;
dist1 = makeAlphaEtaMuFormatI(alpha1,eta1,mu1,Omega1,seriesTol);
dist2 = makeAlphaEtaMuFormatI(alpha2,eta2,mu2,Omega2,seriesTol);

NsSet = 1:NT-1;
NdSet = NT-NsSet;
nAlloc = numel(NsSet);

% Start timing the analytical optimization only. The later Monte Carlo
% validation is deliberately excluded from this runtime.
tOptimization = tic;

%% Unconstrained PSR optimum for each candidate-position allocation
lambdaStarByNs = nan(1,nAlloc);
TstarByNs = nan(1,nAlloc);
TfixedByNs = nan(1,nAlloc);

% A global grid followed by local refinement avoids relying on a
% convexity assumption, without assuming convexity.
lambdaCoarse = linspace(0.05,0.92,49);

fprintf('\nFigure 13: joint optimization versus baselines\n');
fprintf('  NT = %d, Ebar_th = %.1f, gammaBar = %.1f dB\n', ...
    NT,EbarTh,gammaBarDb);
fprintf('  Finding the unconstrained optimum for each allocation ...\n');

for n = 1:nAlloc
    Ns = NsSet(n);
    Nd = NdSet(n);
    Tcoarse = zeros(size(lambdaCoarse));
    for k = 1:numel(lambdaCoarse)
        Tcoarse(k) = exactPsrThroughput(lambdaCoarse(k),R,zeta, ...
            gammaBar,Ns,Nd,dist1,dist2);
    end
    [lambdaStarByNs(n),TstarByNs(n)] = refineMaximum( ...
        lambdaCoarse,Tcoarse,R,zeta,gammaBar,Ns,Nd,dist1,dist2);
    TfixedByNs(n) = exactPsrThroughput(lambdaFixed,R,zeta, ...
        gammaBar,Ns,Nd,dist1,dist2);

    fprintf(['    (Ns,Nd)=(%d,%d): lambda*=%.6f, ', ...
             'T*=%.6f, T(lambda=0.6)=%.6f\n'], ...
        Ns,Nd,lambdaStarByNs(n),TstarByNs(n),TfixedByNs(n));
end

%% Exact feasibility limits for the three baselines
% For a fixed lambda, P_E = F_Z1(x_E)^Ns with
% x_E = 2*Ebar_th/(zeta*lambda*gammaBar).
xEnergyFixed = 2*EbarTh/(zeta*lambdaFixed*gammaBar);
parentProbFixed = parentCdf(xEnergyFixed,dist1);

epsilonNeitherLimit = parentProbFixed^4;
epsilonAllocationLimit = parentProbFixed^(NT-1);

% Protocol-only uses Ns=4 and remains feasible while lambda_min<1.
% Its limiting value is approached as lambda tends to one.
xEnergyLambdaOne = 2*EbarTh/(zeta*gammaBar);
epsilonProtocolLimit = parentCdf(xEnergyLambdaOne,dist1)^4;

fprintf('\n  Reliability feasibility limits\n');
fprintf('    Neither baseline:        epsilon_E = %.6e\n', ...
    epsilonNeitherLimit);
fprintf('    Allocation-only:         epsilon_E = %.6e\n', ...
    epsilonAllocationLimit);
fprintf('    Protocol-only (limit):   epsilon_E = %.6e\n', ...
    epsilonProtocolLimit);

%% Optimize the four strategies over epsilon_E
% Include the exact feasibility limits and a point immediately above the
% open protocol-only limit so that every curve terminates correctly.
epsilonNearProtocolLimit = epsilonProtocolLimit*(1+1e-7);
epsilonGrid = unique([logspace(-6,-1,151),epsilonNeitherLimit, ...
    epsilonAllocationLimit,epsilonProtocolLimit, ...
    epsilonNearProtocolLimit]);
epsilonGrid = epsilonGrid(epsilonGrid >= 1e-6 & epsilonGrid <= 1e-1);

nEpsilon = numel(epsilonGrid);
Tjoint = nan(1,nEpsilon);
Tallocation = nan(1,nEpsilon);
Tprotocol = nan(1,nEpsilon);
Tneither = nan(1,nEpsilon);

NsJoint = nan(1,nEpsilon);
NdJoint = nan(1,nEpsilon);
lambdaJoint = nan(1,nEpsilon);
NsAllocation = nan(1,nEpsilon);

lambdaTolerance = 1e-11;

fprintf('\n  Evaluating the reliability sweep ...\n');
for e = 1:nEpsilon
    epsilonE = epsilonGrid(e);

    % Joint optimization over all feasible allocations and lambda.
    bestT = -Inf;
    for n = 1:nAlloc
        Ns = NsSet(n);
        Nd = NdSet(n);
        qNs = parentQuantile(epsilonE^(1/Ns),dist1);
        lambdaMin = 2*EbarTh/(zeta*gammaBar*qNs);

        if lambdaMin < 1-lambdaTolerance
            if lambdaMin <= lambdaStarByNs(n)
                lambdaCandidate = lambdaStarByNs(n);
                TCandidate = TstarByNs(n);
            else
                lambdaCandidate = lambdaMin;
                TCandidate = exactPsrThroughput(lambdaCandidate,R,zeta, ...
                    gammaBar,Ns,Nd,dist1,dist2);
            end

            if TCandidate > bestT
                bestT = TCandidate;
                Tjoint(e) = TCandidate;
                NsJoint(e) = Ns;
                NdJoint(e) = Nd;
                lambdaJoint(e) = lambdaCandidate;
            end
        end
    end

    % Allocation-only: lambda is fixed and Ns is enumerated.
    bestT = -Inf;
    for n = 1:nAlloc
        Ns = NsSet(n);
        qNs = parentQuantile(epsilonE^(1/Ns),dist1);
        lambdaMin = 2*EbarTh/(zeta*gammaBar*qNs);
        if lambdaMin <= lambdaFixed+lambdaTolerance && ...
                TfixedByNs(n) > bestT
            bestT = TfixedByNs(n);
            Tallocation(e) = TfixedByNs(n);
            NsAllocation(e) = Ns;
        end
    end

    % Protocol-only: the balanced allocation is fixed at (4,4).
    nBalanced = 4;
    qBalanced = parentQuantile(epsilonE^(1/nBalanced),dist1);
    lambdaMinBalanced = 2*EbarTh/(zeta*gammaBar*qBalanced);
    if lambdaMinBalanced < 1-lambdaTolerance
        if lambdaMinBalanced <= lambdaStarByNs(nBalanced)
            lambdaProtocol = lambdaStarByNs(nBalanced);
            Tprotocol(e) = TstarByNs(nBalanced);
        else
            lambdaProtocol = lambdaMinBalanced;
            Tprotocol(e) = exactPsrThroughput(lambdaProtocol,R,zeta, ...
                gammaBar,4,4,dist1,dist2);
        end
    end

    % Neither variable is optimized: (Ns,Nd)=(4,4), lambda=0.6.
    if lambdaMinBalanced <= lambdaFixed+lambdaTolerance
        Tneither(e) = TfixedByNs(nBalanced);
    end
end

% Stop the analytical-optimization timer BEFORE Monte Carlo validation.
optimizationRuntime = toc(tOptimization);

fprintf('\nAnalytical optimization runtime = %.3f s\n', ...
    optimizationRuntime);

% Save a compact runtime record that can be uploaded/shared directly.
runtimeFile = 'Fig13_Runtime_Info.txt';
fidRuntime = fopen(runtimeFile,'w');
if fidRuntime ~= -1
    fprintf(fidRuntime,'Figure 13 analytical optimization runtime\n');
    fprintf(fidRuntime,'MATLAB version: %s\n',matlabVersionText);
    fprintf(fidRuntime,'MATLAB release: %s\n',matlabReleaseText);
    fprintf(fidRuntime,'Architecture: %s\n',architectureText);
    fprintf(fidRuntime,'CPU: %s\n',cpuText);
    fprintf(fidRuntime,'NT: %d\n',NT);
    fprintf(fidRuntime,'Coarse PSR grid points: %d\n',numel(lambdaCoarse));
    fprintf(fidRuntime,'Coarse PSR grid interval: [%.2f, %.2f]\n', ...
        lambdaCoarse(1),lambdaCoarse(end));
    fprintf(fidRuntime,'Local-search TolX: 1e-8\n');
    fprintf(fidRuntime,'Feasibility tolerance: %.1e\n',lambdaTolerance);
    fprintf(fidRuntime,'Inverse-CDF bisection iterations: 90\n');
    fprintf(fidRuntime,'Analytical optimization runtime: %.6f s\n', ...
        optimizationRuntime);
    fprintf(fidRuntime,'Monte Carlo validation excluded from runtime: yes\n');
    fclose(fidRuntime);
    fprintf('Runtime summary saved to %s\n',runtimeFile);
else
    warning('Could not create %s.',runtimeFile);
end

%% Key numerical values used in the discussion
epsilonReport = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6];
reportIndex = nearestLogIndices(epsilonGrid,epsilonReport);

fprintf('\n  Key optimized results\n');
for k = 1:numel(epsilonReport)
    e = reportIndex(k);
    fprintf(['    epsilon_E=%.0e: joint T=%.6f, Ns*=%g, ', ...
             'lambda*=%.6f; allocation-only T=%.6f; ', ...
             'protocol-only T=%.6f; neither T=%.6f\n'], ...
        epsilonReport(k),Tjoint(e),NsJoint(e),lambdaJoint(e), ...
        Tallocation(e),Tprotocol(e),Tneither(e));
end

%% Monte Carlo validation at one point per decade
epsilonMC = epsilonReport;
mcIndex = reportIndex;

NsDesignMC = nan(4,numel(epsilonMC));
NdDesignMC = nan(4,numel(epsilonMC));
lambdaDesignMC = nan(4,numel(epsilonMC));

% Row 1: joint design.
NsDesignMC(1,:) = NsJoint(mcIndex);
NdDesignMC(1,:) = NdJoint(mcIndex);
lambdaDesignMC(1,:) = lambdaJoint(mcIndex);

% Row 2: allocation-only design.
NsDesignMC(2,:) = NsAllocation(mcIndex);
NdDesignMC(2,:) = NT-NsAllocation(mcIndex);
lambdaDesignMC(2,isfinite(NsAllocation(mcIndex))) = lambdaFixed;

% Row 3: protocol-only design.
protocolFeasibleMC = isfinite(Tprotocol(mcIndex));
NsDesignMC(3,protocolFeasibleMC) = 4;
NdDesignMC(3,protocolFeasibleMC) = 4;
for k = find(protocolFeasibleMC)
    epsilonE = epsilonMC(k);
    qBalanced = parentQuantile(epsilonE^(1/4),dist1);
    lambdaMinBalanced = 2*EbarTh/(zeta*gammaBar*qBalanced);
    lambdaDesignMC(3,k) = max(lambdaStarByNs(4),lambdaMinBalanced);
end

% Row 4: neither variable is optimized.
neitherFeasibleMC = isfinite(Tneither(mcIndex));
NsDesignMC(4,neitherFeasibleMC) = 4;
NdDesignMC(4,neitherFeasibleMC) = 4;
lambdaDesignMC(4,neitherFeasibleMC) = lambdaFixed;

nMC = 1e6;
chunkSize = 1e5;
fprintf('\n  Running Monte Carlo validation with %d realizations ...\n',nMC);
Tsimulation = monteCarloDesignThroughput(NsDesignMC,NdDesignMC, ...
    lambdaDesignMC,nMC,chunkSize,R,zeta,gammaBar,NT, ...
    alpha1,eta1,mu1,Omega1,alpha2,eta2,mu2,Omega2);

TanalyticalAtMC = [Tjoint(mcIndex);Tallocation(mcIndex); ...
                   Tprotocol(mcIndex);Tneither(mcIndex)];
for strategy = 1:4
    fprintf('    Maximum |MC-exact|, strategy %d = %.3e\n', ...
        strategy,maxFinite(abs(Tsimulation(strategy,:)- ...
        TanalyticalAtMC(strategy,:))));
end

%% Plot: single-column layout consistent with Figs. 1-10
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

blue = [0.0000,0.4470,0.7410];
orange = [0.8500,0.3250,0.0980];
green = [0.4660,0.6740,0.1880];
purple = [0.4940,0.1840,0.5560];

hJoint = plot(ax,epsilonGrid,Tjoint,'-','Color',blue, ...
    'LineWidth',1.2);
hAllocation = plot(ax,epsilonGrid,Tallocation,'--','Color',orange, ...
    'LineWidth',1.2);
hProtocol = plot(ax,epsilonGrid,Tprotocol,'-.','Color',green, ...
    'LineWidth',1.2);
hNeither = plot(ax,epsilonGrid,Tneither,':','Color',purple, ...
    'LineWidth',1.2);

% Analytical curves are drawn first.  Transparent open black circles then
% show the simulations without hiding the curves underneath.
hMC = plot(ax,epsilonMC,Tsimulation(1,:),'ko','LineStyle','none', ...
    'MarkerFaceColor','none','MarkerSize',5,'LineWidth',0.8);
for strategy = 2:4
    plot(ax,epsilonMC,Tsimulation(strategy,:),'ko', ...
        'LineStyle','none','MarkerFaceColor','none','MarkerSize',5, ...
        'LineWidth',0.8,'HandleVisibility','off');
end

% Identical black crosses mark the reliability limits beyond which a
% baseline is infeasible. Infeasibility is never plotted as zero.
hFeasibility = plot(ax,epsilonAllocationLimit,TfixedByNs(7),'kx', ...
    'MarkerSize',7,'LineWidth',1.0,'LineStyle','none');
plot(ax,epsilonProtocolLimit,0,'kx', ...
    'MarkerSize',7,'LineWidth',1.0,'LineStyle','none', ...
    'HandleVisibility','off');
plot(ax,epsilonNeitherLimit,TfixedByNs(4),'kx', ...
    'MarkerSize',7,'LineWidth',1.0,'LineStyle','none', ...
    'HandleVisibility','off');

xlabel(ax,'Energy-outage target, $\epsilon_E$', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');
ylabel(ax,'Optimised throughput, $\mathcal{T}_{\mathrm{PSR}}$ (bit/s/Hz)', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

set(ax, ...
    'XScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',0.8, ...
    'TickDir','in', ...
    'XMinorGrid','off', ...
    'YMinorGrid','on', ...
    'GridAlpha',0.22, ...
    'MinorGridAlpha',0.12, ...
    'Layer','top');
xlim(ax,[1e-6,1e-1]);
ylim(ax,[0,2.08]);
set(ax,'XTick',10.^(-6:-1));

% A compact two-column legend leaves room for the magnified view.
leg = legend(ax, ...
    [hJoint,hAllocation,hProtocol,hNeither,hMC,hFeasibility], ...
    {'Joint','Allocation-only','Protocol-only','Fully fixed', ...
     'Simulation','Feasibility limit'}, ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'NumColumns',2, ...
    'Box','on', ...
    'Color','w', ...
    'EdgeColor','k', ...
    'Interpreter','tex');

leg.LineWidth = 1;
leg.ItemTokenSize = [10 10];
set(leg,'Units','normalized','Position',[0.39,0.12,0.54,0.22]);

%% Magnified view of the overlapping high-throughput region
% This inset uses the same quantities and curves as the main axes.  It is
% therefore a true magnification, not an additional design-variable plot.
% Simulation markers and a second legend are omitted to keep it readable.
axZoom = axes(fig,'Position',[0.4,0.52,0.47,0.30]);
hold(axZoom,'on');
box(axZoom,'on');
grid(axZoom,'on');

plot(axZoom,epsilonGrid,Tjoint,'-','Color',blue, ...
    'LineWidth',1.0);
plot(axZoom,epsilonGrid,Tallocation,'--','Color',orange, ...
    'LineWidth',1.0);
plot(axZoom,epsilonGrid,Tprotocol,'-.','Color',green, ...
    'LineWidth',1.0);
plot(axZoom,epsilonGrid,Tneither,':','Color',purple, ...
    'LineWidth',1.0);

set(axZoom, ...
    'XScale','log', ...
    'XLim',[1e-5,1e-3], ...
    'YLim',[1.90,2.01], ...
    'XTick',[1e-5,1e-4,1e-3], ...
    'YTick',[1.90,1.95,2.00], ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'LineWidth',0.65, ...
    'TickDir','in', ...
    'XMinorGrid','off', ...
    'YMinorGrid','on', ...
    'GridAlpha',0.22, ...
    'MinorGridAlpha',0.12, ...
    'Layer','top');

%% Export publication files and numerical tables
outputDir = fullfile(pwd,'Fig13_output');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

if exist('exportgraphics','file') == 2
    exportgraphics(fig,fullfile(outputDir, ...
        'Fig13_joint_optimization.pdf'),'ContentType','vector');
    exportgraphics(fig,fullfile(outputDir, ...
        'Fig13_joint_optimization.png'),'Resolution',600);
else
    print(fig,fullfile(outputDir, ...
        'Fig13_joint_optimization.pdf'),'-dpdf','-painters');
    print(fig,fullfile(outputDir, ...
        'Fig13_joint_optimization.png'),'-dpng','-r600');
end
savefig(fig,fullfile(outputDir,'Fig13_joint_optimization.fig'));

curveTable = table(epsilonGrid(:),Tjoint(:),Tallocation(:), ...
    Tprotocol(:),Tneither(:),NsJoint(:),NdJoint(:), ...
    lambdaJoint(:),NsAllocation(:), ...
    'VariableNames',{'epsilon_E','T_joint','T_allocation_only', ...
    'T_protocol_only','T_neither','Ns_joint','Nd_joint', ...
    'lambda_joint','Ns_allocation_only'});
writetable(curveTable,fullfile(outputDir,'Fig13_curves.csv'));

mcTable = table(epsilonMC(:),Tsimulation(1,:).', ...
    Tsimulation(2,:).',Tsimulation(3,:).',Tsimulation(4,:).', ...
    'VariableNames',{'epsilon_E','T_joint_MC','T_allocation_MC', ...
    'T_protocol_MC','T_neither_MC'});
writetable(mcTable,fullfile(outputDir,'Fig13_monte_carlo.csv'));

allocationTable = table(NsSet(:),NdSet(:),lambdaStarByNs(:), ...
    TstarByNs(:),TfixedByNs(:), ...
    'VariableNames',{'Ns','Nd','lambda_unconstrained', ...
    'T_unconstrained','T_at_lambda_0p6'});
writetable(allocationTable,fullfile(outputDir, ...
    'Fig13_allocation_summary.csv'));

%% Local functions
function value = maxFinite(x)
% Maximum after excluding NaN and Inf values.
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = max(x);
    end
end

function indices = nearestLogIndices(grid,values)
% Indices of grid points nearest to specified positive values.
    indices = zeros(size(values));
    for k = 1:numel(values)
        [~,indices(k)] = min(abs(log10(grid)-log10(values(k))));
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
% Parent alpha-eta-mu channel-power CDF.
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

function T = exactPsrThroughput(lambda,R,zeta,gammaBar, ...
    Ns,Nd,dist1,dist2)
% Exact PSR delay-limited throughput from the one-dimensional outage representation.
    if ~isfinite(lambda) || lambda <= 0 || lambda >= 1
        T = 0;
        return;
    end

    a = (1-lambda)*gammaBar;
    b = zeta*lambda*gammaBar;
    gammaTh = 2^(2*R)-1;
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

function [lambdaStar,TStar] = refineMaximum(grid,Tgrid,R,zeta, ...
    gammaBar,Ns,Nd,dist1,dist2)
% Refine the best coarse-grid point using a bounded local search.
    [~,index] = max(Tgrid);
    indexLower = max(1,index-1);
    indexUpper = min(numel(grid),index+1);
    lower = max(grid(indexLower),1e-6);
    upper = min(grid(indexUpper),1-1e-6);
    objective = @(lambda) -exactPsrThroughput(lambda,R,zeta, ...
        gammaBar,Ns,Nd,dist1,dist2);
    options = optimset('TolX',1e-8,'Display','off');
    lambdaStar = fminbnd(objective,lower,upper,options);
    TStar = -objective(lambdaStar);
end

function Tsim = monteCarloDesignThroughput(NsDesign,NdDesign, ...
    lambdaDesign,nMC,chunkSize,R,zeta,gammaBar,NT, ...
    alpha1,eta1,mu1,Omega1,alpha2,eta2,mu2,Omega2)
% Monte Carlo simulation of the end-to-end AF SNR for all displayed designs.
    [nStrategy,nPoint] = size(lambdaDesign);
    success = zeros(nStrategy,nPoint);
    completed = 0;
    gammaTh = 2^(2*R)-1;

    while completed < nMC
        nNow = min(chunkSize,nMC-completed);
        Z1 = sampleAlphaEtaMuFormatI( ...
            nNow,NT-1,alpha1,eta1,mu1,Omega1);
        Z2 = sampleAlphaEtaMuFormatI( ...
            nNow,NT-1,alpha2,eta2,mu2,Omega2);

        for strategy = 1:nStrategy
            for point = 1:nPoint
                lambda = lambdaDesign(strategy,point);
                Ns = NsDesign(strategy,point);
                Nd = NdDesign(strategy,point);
                if isfinite(lambda) && isfinite(Ns) && isfinite(Nd)
                    X = max(Z1(:,1:Ns),[],2);
                    Y = max(Z2(:,1:Nd),[],2);
                    a = (1-lambda)*gammaBar;
                    b = zeta*lambda*gammaBar;
                    gammaAF = (a*b*X.^2.*Y)./(a*X+b*X.*Y+1);
                    success(strategy,point) = ...
                        success(strategy,point)+sum(gammaAF >= gammaTh);
                end
            end
        end
        completed = completed+nNow;
    end

    Tsim = R*success/nMC;
    Tsim(~isfinite(lambdaDesign)) = NaN;
end

function Z = sampleAlphaEtaMuFormatI(nRows,nCols,alpha,eta,mu,Omega)
% Exact Format-I sampler for arbitrary real mu > 0.
    sigmaY2 = Omega^(alpha/2)/(2*mu*(1+eta));
    sigmaX2 = eta*sigmaY2;
    Gx = (2*sigmaX2)*randg(mu,nRows,nCols);
    Gy = (2*sigmaY2)*randg(mu,nRows,nCols);
    Z = (Gx+Gy).^(2/alpha);
end
