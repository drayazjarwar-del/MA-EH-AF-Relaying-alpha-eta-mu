%% Figure 11: Energy-reliability feasibility conditions
% Implements the energy-reliability feasibility conditions for PSR over Format-I alpha-eta-mu fading.
%
% The figure plots lambda_min(Ns). Separate TSR curves are unnecessary
% because rho_min(Ns) = lambda_min(Ns)/2 exactly. Moreover, the horizontal
% lambda_min = 1 boundary also provides the IRR feasibility test because
% chi_IRR equals chi_PSR evaluated at lambda = 1.

clear; close all; clc;
rng(9,'twister');
fprintf('\nRunning Fig11_Energy_Feasibility.m (analytical + Monte Carlo)\n');

%% Common fading and system parameters
alpha = 2;
eta   = 0.5;       % Format I
mu    = 2;
Omega = 1;

zeta       = 0.8;
snrDb      = 20;
gammaBar   = 10^(snrDb/10);  % P_S/N_0
EbarTh     = 20;              % E_th/(N_0*T)
epsilonE   = [1e-1, 1e-2, 1e-3];
NsValues   = 1:12;

%% Numerical controls
seriesTol = 1e-12;
inverseTol = 1e-12;
nMC = 1e6;
chunkSize = 1e5;

% Construct the positive gamma-mixture representation of the parent CDF once.
cdfModel = makeAlphaEtaMuCdfModel( ...
    alpha, eta, mu, Omega, seriesTol);

%% Quantiles and minimum feasible protocol factors
nTarget = numel(epsilonE);
nNs = numel(NsValues);
qNs = zeros(nTarget,nNs);
lambdaMin = zeros(nTarget,nNs);
rhoMin = zeros(nTarget,nNs);
NsMin = NaN(1,nTarget);

for e = 1:nTarget
    for j = 1:nNs
        Ns = NsValues(j);
        probabilityLevel = epsilonE(e)^(1/Ns);

        qNs(e,j) = inverseAlphaEtaMuCdf( ...
            probabilityLevel, cdfModel, inverseTol);

        lambdaMin(e,j) = 2*EbarTh/(zeta*gammaBar*qNs(e,j));
        rhoMin(e,j) = lambdaMin(e,j)/2;
    end

    firstFeasible = find(lambdaMin(e,:) < 1,1,'first');
    if ~isempty(firstFeasible)
        NsMin(e) = NsValues(firstFeasible);
    end
end

assert(isequal(NsMin,[2,3,4]), ...
    'Expected N_s,min = [2,3,4] for the selected reliability targets.');

%% Monte Carlo feasibility boundaries
% For each Ns, simulate the maximum selected first-hop gain. Its empirical
% epsilon_E-quantile is the Monte Carlo counterpart of the analytical q_Ns quantile.
% Substitution into the feasibility expression gives the simulated lambda_min marker.
[lambdaMinMC,qSelectedMC,outageAtExactBoundary] = ...
    monteCarloLambdaMin(nMC,chunkSize,NsValues,epsilonE, ...
    alpha,eta,mu,Omega,EbarTh,zeta,gammaBar,lambdaMin);

relativeMCError = abs(lambdaMinMC-lambdaMin)./lambdaMin;
targetMatrix = repmat(epsilonE(:),1,nNs);
boundaryOutageError = abs(outageAtExactBoundary-targetMatrix);

%% Display and save numerical values
resultsTable = table(NsValues(:), ...
    lambdaMin(1,:).',lambdaMin(2,:).',lambdaMin(3,:).', ...
    lambdaMinMC(1,:).',lambdaMinMC(2,:).',lambdaMinMC(3,:).', ...
    rhoMin(1,:).',rhoMin(2,:).',rhoMin(3,:).', ...
    outageAtExactBoundary(1,:).',outageAtExactBoundary(2,:).', ...
    outageAtExactBoundary(3,:).', ...
    'VariableNames',{'Ns', ...
    'lambdaMin_eps1e_1','lambdaMin_eps1e_2','lambdaMin_eps1e_3', ...
    'lambdaMinMC_eps1e_1','lambdaMinMC_eps1e_2','lambdaMinMC_eps1e_3', ...
    'rhoMin_eps1e_1','rhoMin_eps1e_2','rhoMin_eps1e_3', ...
    'MCoutageAtExact_eps1e_1','MCoutageAtExact_eps1e_2', ...
    'MCoutageAtExact_eps1e_3'});

disp(resultsTable);
fprintf('Figure 9 diagnostics\n');
fprintf('  alpha = %.2f, eta = %.2f, mu = %.2f, Omega = %.2f\n', ...
    alpha,eta,mu,Omega);
fprintf('  gammaBar = %.0f dB, EbarTh = %.0f\n',snrDb,EbarTh);
fprintf('  CDF truncation order K = %d\n',cdfModel.K);
fprintf('  Remaining mixture-weight residual = %.3e\n', ...
    cdfModel.residual);
fprintf('  N_s,min = [%d, %d, %d]\n',NsMin);
fprintf('  Monte Carlo realizations = %d\n',nMC);
fprintf('  Maximum relative difference: MC versus analytical lambda_min = %.3e\n', ...
    max(relativeMCError(:)));
fprintf('  Maximum absolute energy-outage error at analytical boundary = %.3e\n', ...
    max(boundaryOutageError(:)));

%% Plot
% Publication style used consistently in Figures 2-8.
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

yLower = 0.30;
yUpper = 5.50;

% The region lambda_min >= 1 is infeasible for PSR because lambda < 1.
hShade = patch(ax, ...
    [NsValues(1)-0.45, NsValues(end)+0.45, ...
     NsValues(end)+0.45, NsValues(1)-0.45], ...
    [1,1,yUpper,yUpper], ...
    [0.82,0.82,0.82], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.40, ...
    'HandleVisibility','off');

% Same blue-orange-green palette used in the preceding figures.
colors = [
    0.0000 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.4660 0.6740 0.1880
];
styles = {'-','--','-.'};
hCurve = gobjects(1,nTarget);

for e = 1:nTarget
    hCurve(e) = semilogy(ax,NsValues,lambdaMin(e,:), ...
        'LineStyle',styles{e},'Color',colors(e,:), ...
        'LineWidth',1.2);

    % Match Figures 2-8: unfilled black circles denote simulation.
    semilogy(ax,NsValues,lambdaMinMC(e,:),'ko', ...
        'LineStyle','none','MarkerSize',5,'LineWidth',0.8, ...
        'MarkerFaceColor','none','HandleVisibility','off');

    jMin = find(NsValues==NsMin(e),1);
    semilogy(ax,NsMin(e),lambdaMin(e,jMin),'p', ...
    'Color','k', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',9, ...
    'LineWidth',0.8, ...
    'HandleVisibility','off');
end

hBoundary = yline(ax,1,'k--','LineWidth',1.0);
hSimulation = semilogy(ax,NaN,NaN,'ko', ...
    'MarkerFaceColor','none','MarkerSize',5, ...
    'LineWidth',0.6,'LineStyle','none');
hMinMarker = semilogy(ax,NaN,NaN,'kp', ...
    'MarkerFaceColor','k','MarkerSize',8, ...
    'LineWidth',0.6,'LineStyle','none');

text(ax,6.5,1.2,'PSR infeasible', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'HorizontalAlignment','center');

xlabel(ax,'Number of source candidate positions, $N_s$', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

ylabel(ax,'Minimum power-splitting factor, $\lambda_{\min}(N_s)$', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

xlim(ax,[NsValues(1)-0.45,NsValues(end)+0.45]);
ylim(ax,[yLower,yUpper]);
xticks(ax,NsValues);
yticks(ax,[0.3,0.4,0.5,0.7,1,2,3,5]);
yticklabels(ax,{'0.3','0.4','0.5','0.7','1','2','3','5'});

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

lgd = legend(ax,[hCurve,hSimulation,hMinMarker,hBoundary], ...
    {'$\epsilon_E=10^{-1}\;(N_{s,\min}=2)$', ...
     '$\epsilon_E=10^{-2}\;(N_{s,\min}=3)$', ...
     '$\epsilon_E=10^{-3}\;(N_{s,\min}=4)$', ...
     'Simulation','First feasible $N_s$','PSR feasibility boundary'}, ...
    'Interpreter','latex', ...
    'Location','northeast', ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'Box','on');

lgd.LineWidth = 1;
lgd.ItemTokenSize = [11 11];

% Add only a small amount of width towards the left while keeping the
% right edge fixed.  This prevents clipping without leaving excess space.
drawnow;
lgd.Units = 'normalized';
legendPosition = lgd.Position;
legendRightEdge = legendPosition(1)+legendPosition(3);
legendPosition(1) = max(0.02,legendPosition(1)-0.04);
legendPosition(3) = legendRightEdge-legendPosition(1);
lgd.Position = legendPosition;

%% Export publication files and numerical table
outputDir = fullfile(pwd,'Fig11_output');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

writetable(resultsTable, ...
    fullfile(outputDir,'Fig11_Energy_Feasibility_values.csv'));

if exist('exportgraphics','file') == 2
    exportgraphics(fig, ...
        fullfile(outputDir,'Fig11_Energy_Feasibility.pdf'), ...
        'ContentType','vector');
    exportgraphics(fig, ...
        fullfile(outputDir,'Fig11_Energy_Feasibility.png'), ...
        'Resolution',600);
else
    print(fig,fullfile(outputDir,'Fig11_Energy_Feasibility.pdf'), ...
        '-dpdf','-painters');
    print(fig,fullfile(outputDir,'Fig11_Energy_Feasibility.png'), ...
        '-dpng','-r600');
end
savefig(fig,fullfile(outputDir,'Fig11_Energy_Feasibility.fig'));

%% Local functions
function model = makeAlphaEtaMuCdfModel(alpha,eta,mu,Omega,tol)
% Positive gamma-mixture representation of the Format-I CDF in Eq. (18).
% MATLAB's regularized gammainc is used, so Gamma(shape) is absorbed into
% each normalized mixture weight.

    validateattributes(alpha,{'numeric'},{'scalar','real','positive'});
    validateattributes(eta,  {'numeric'},{'scalar','real','positive'});
    validateattributes(mu,   {'numeric'},{'scalar','real','positive'});
    validateattributes(Omega,{'numeric'},{'scalar','real','positive'});

    h = (1+eta)^2/(4*eta);
    H = (1-eta^2)/(4*eta);

    maxK = 100000;
    weights = zeros(1,64);
    weights(1) = h^(-mu);
    sumWeights = weights(1);
    K = 0;

    while (1-sumWeights) > tol
        if K >= maxK
            error('makeAlphaEtaMuCdfModel:NoConvergence', ...
                'The CDF series did not meet the requested tolerance.');
        end
        if K+2 > numel(weights)
            weights = [weights,zeros(1,numel(weights))]; %#ok<AGROW>
        end
        ratio = (H/h)^2*(mu+K)/(K+1);
        weights(K+2) = weights(K+1)*ratio;
        K = K+1;
        sumWeights = sumWeights+weights(K+1);
    end

    model.alpha = alpha;
    model.mu = mu;
    model.Omega = Omega;
    model.h = h;
    model.weights = weights(1:K+1);
    model.shapes = 2*(mu+(0:K));
    model.K = K;
    model.residual = max(0,1-sumWeights);
end

function F = alphaEtaMuCdfFormatI(z,model)
% Evaluate the parent power-gain CDF using the precomputed model.

    z = max(z,0);
    x = 2*model.mu*model.h*(z./model.Omega).^(model.alpha/2);
    F = zeros(size(z));
    for k = 1:numel(model.weights)
        F = F+model.weights(k)*gammainc(x,model.shapes(k),'lower');
    end
    F = min(max(F,0),1);
end

function q = inverseAlphaEtaMuCdf(p,model,tol)
% Generalized inverse F^{-1}(p), evaluated by bracket expansion followed
% by bisection. The parent CDF is continuous and strictly increasing.

    validateattributes(p,{'numeric'},{'scalar','real','>',0,'<',1});

    lower = 0;
    upper = model.Omega;
    while alphaEtaMuCdfFormatI(upper,model) < p
        upper = 2*upper;
        if ~isfinite(upper)
            error('inverseAlphaEtaMuCdf:NoBracket', ...
                'Unable to bracket the requested CDF probability.');
        end
    end

    maxIterations = 250;
    for iteration = 1:maxIterations
        midpoint = 0.5*(lower+upper);
        if alphaEtaMuCdfFormatI(midpoint,model) < p
            lower = midpoint;
        else
            upper = midpoint;
        end

        if (upper-lower) <= tol*max(1,midpoint)
            break;
        end
    end
    q = 0.5*(lower+upper);
end

function [lambdaMC,qMC,pAtExact] = monteCarloLambdaMin( ...
    nMC,chunkSize,NsValues,epsilonE,alpha,eta,mu,Omega, ...
    EbarTh,zeta,gammaBar,lambdaExact)
% Monte Carlo validation of the feasibility boundary.
%
% For each channel realization and each Ns, Xselected stores the maximum
% of Ns independent candidate-position gains. The empirical epsilon_E
% quantile of Xselected is then inserted into the feasibility expression.

    NsMax = max(NsValues);
    nNs = numel(NsValues);
    nTarget = numel(epsilonE);

    % About 96 MB for nMC=1e6 and NsMax=12. Channel generation is chunked
    % so that the complete parent-gain arrays are never held simultaneously.
    Xselected = zeros(nMC,nNs);

    nCompleted = 0;
    while nCompleted < nMC
        nNow = min(chunkSize,nMC-nCompleted);
        rows = nCompleted+(1:nNow);

        Z = sampleAlphaEtaMuFormatI( ...
            nNow,NsMax,alpha,eta,mu,Omega);
        Zmax = cummax(Z,2);
        Xselected(rows,:) = Zmax(:,NsValues);

        nCompleted = nCompleted+nNow;
    end

    qMC = zeros(nTarget,nNs);
    lambdaMC = zeros(nTarget,nNs);
    pAtExact = zeros(nTarget,nNs);

    for j = 1:nNs
        sortedSelected = sort(Xselected(:,j),'ascend');

        for e = 1:nTarget
            qMC(e,j) = empiricalQuantileSorted( ...
                sortedSelected,epsilonE(e));
            lambdaMC(e,j) = 2*EbarTh/( ...
                zeta*gammaBar*qMC(e,j));

            exactThreshold = 2*EbarTh/( ...
                zeta*gammaBar*lambdaExact(e,j));
            pAtExact(e,j) = mean(Xselected(:,j) < exactThreshold);
        end
    end
end

function Z = sampleAlphaEtaMuFormatI(nRows,nCols,alpha,eta,mu,Omega)
% Exact Format-I alpha-eta-mu power-gain sampler for arbitrary real mu>0.
% The generalized scale satisfies Omega=(E[Z^(alpha/2)])^(2/alpha).

    sigmaY2 = Omega^(alpha/2)/(2*mu*(1+eta));
    sigmaX2 = eta*sigmaY2;

    Gx = (2*sigmaX2)*randg(mu,nRows,nCols);
    Gy = (2*sigmaY2)*randg(mu,nRows,nCols);
    Z = (Gx+Gy).^(2/alpha);
end

function q = empiricalQuantileSorted(sortedSamples,p)
% Linearly interpolated empirical quantile using r=1+(n-1)p.

    n = numel(sortedSamples);
    r = 1+(n-1)*p;
    lowerIndex = floor(r);
    upperIndex = ceil(r);

    if lowerIndex == upperIndex
        q = sortedSamples(lowerIndex);
    else
        fraction = r-lowerIndex;
        q = (1-fraction)*sortedSamples(lowerIndex)+ ...
            fraction*sortedSamples(upperIndex);
    end
end
