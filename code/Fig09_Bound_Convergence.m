%% Figure 9: Convergence of finite-sum bounds with finite-partition error control
% Revised version using probability-mass midpoint refinement
%
% This script generates Fig. 9 for the PSR protocol using the original
% numerical parameters of the paper:
%
%   alpha_1 = alpha_2 = 2
%   eta_1   = eta_2   = 0.5
%   mu_1    = mu_2    = 2
%   Ns = Nd = 3
%   zeta = 0.8
%   lambda = 0.6
%   R = 2 bit/s/Hz
%
% The plotted metric is
%
%   epsilon_L =
%      (P_UB,L - P_LB,L) /
%      [0.5*(P_UB,L + P_LB,L)].
%
% The adaptive algorithm repeatedly selects the interval with the largest
% local finite-partition gap and divides it at the midpoint of its first-hop
% probability mass:
%
%   F_X(x_mid) = [F_X(x_left)+F_X(x_right)]/2.
%
% This refinement guarantees that the finite-partition gap contribution of the
% selected interval is reduced by one half.
%
% Exported files:
%   Fig09_Bound_Convergence.pdf
%   Fig09_Bound_Convergence.png
%   Fig09_Bound_Convergence_Data.csv
%   Fig09_Bound_Convergence_Data.mat
%
% Public release filename:
%   Fig09_Bound_Convergence.m

clear;
clc;
close all;

%% ======================================================================
% 1. PAPER PARAMETERS
% =======================================================================

% First hop: S -> R
alpha1  = 2;
eta1    = 0.5;
mu1     = 2;
Omega1  = 1;
format1 = 1;               % 1 = Format I, 2 = Format II
Ns      = 3;

% Second hop: R -> D
alpha2  = 2;
eta2    = 0.5;
mu2     = 2;
Omega2  = 1;
format2 = 1;
Nd      = 3;

% PSR parameters
zeta   = 0.8;
lambda = 0.6;
R      = 2;               % bit/s/Hz

% PSR threshold and protocol constants
gammaTh = 2^(2*R)-1;
A = 1-lambda;
B = zeta*lambda;

% SNR values used in Fig. 7
snrDbCases  = [20 25 30];
snrLinCases = 10.^(snrDbCases/10);

% Alpha-eta-mu mixture settings
seriesTolerance = 1e-14;
maximumTerms    = 500;

% Adaptive convergence settings
targetRelativeGap = 5e-4;
displayTolerance  = 1e-3;
maximumIntervals  = 50000;

%% ========================================================================
% 2. BUILD THE TWO PARENT DISTRIBUTIONS
% ========================================================================

dist1 = buildAlphaEtaMuDistribution( ...
    alpha1,eta1,mu1,Omega1,format1, ...
    seriesTolerance,maximumTerms);

dist2 = buildAlphaEtaMuDistribution( ...
    alpha2,eta2,mu2,Omega2,format2, ...
    seriesTolerance,maximumTerms);

%% ========================================================================
% 3. COMPUTE COMPLETE CONVERGENCE HISTORIES
% ========================================================================

numberOfCases = numel(snrDbCases);
history = cell(numberOfCases,1);

fprintf('Computing revised adaptive-bound convergence histories...\n\n');

for caseIndex = 1:numberOfCases

    gammaBar = snrLinCases(caseIndex);

    a = A*gammaBar;
    b = B*gammaBar;

    history{caseIndex} = adaptiveBoundHistoryProbabilityMidpoint( ...
        a,b,gammaTh,Ns,Nd,dist1,dist2, ...
        targetRelativeGap,maximumIntervals);

    fprintf(['SNR = %2.0f dB | final relative gap = %.3e | ' ...
             'intervals = %d\n'], ...
             snrDbCases(caseIndex), ...
             history{caseIndex}.relativeGap(end), ...
             history{caseIndex}.numberOfIntervals(end));
end

%% ========================================================================
% 4. PUBLICATION-READY FIGURE
% ========================================================================

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

% Same publication style as Figures 2-6
curveColors = [
    0.0000 0.4470 0.7410;   % Blue
    0.8500 0.3250 0.0980;   % Orange
    0.4660 0.6740 0.1880    % Green
];

lineStyles = {'-','--','-.'};

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 3.5 3.10], ...
    'PaperPositionMode','auto');

ax = axes(fig);
hold(ax,'on');
box(ax,'on');
grid(ax,'on');

legendText = cell(numberOfCases,1);
hCurve = gobjects(numberOfCases,1);

for caseIndex = 1:numberOfCases

    hCurve(caseIndex) = loglog( ...
        ax, ...
        history{caseIndex}.numberOfIntervals, ...
        history{caseIndex}.relativeGap, ...
        'LineStyle',lineStyles{caseIndex}, ...
        'Color',curveColors(caseIndex,:), ...
        'LineWidth',1.2);

    legendText{caseIndex} = sprintf( ...
        '%d dB',snrDbCases(caseIndex));
end

% Displayed tolerance line. It is included in the legend so that the
% annotation style is consistent with the other numerical figures.
hTolerance = yline( ...
    ax, ...
    displayTolerance, ...
    ':k', ...
    'LineWidth',1.0, ...
    'DisplayName','$\epsilon_{\mathrm{rel}}=10^{-3}$');

xlabel(ax,'Number of adaptive partition intervals', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

ylabel(ax, ...
    'Relative bound gap, $\Delta_L/\widehat{P}_{\mathrm{out}}^{\,L}$', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

maximumPlottedIntervals = max(cellfun( ...
    @(item) item.numberOfIntervals(end),history));

xlim(ax,[1 maximumPlottedIntervals]);

xticks(ax,[1 10 100 1000 10000]);
xticklabels(ax,{ ...
    '$10^{0}$', ...
    '$10^{1}$', ...
    '$10^{2}$', ...
    '$10^{3}$', ...
    '$10^{4}$'});

set(ax,'TickLabelInterpreter','latex');

allRelativeGaps = cell2mat(cellfun( ...
    @(item) item.relativeGap(:),history,'UniformOutput',false));

minimumGap = min(allRelativeGaps);
minimumExponent = floor(log10(minimumGap))-1;
minimumExponent = min(minimumExponent,-4);

% Place the lower axis slightly below the final convergence points.
finalRelativeGaps = cellfun( ...
    @(item) item.relativeGap(end),history);

lowerAxisLimit = 0.8*min(finalRelativeGaps);

ylim(ax,[lowerAxisLimit 2]);
yticks(ax,[1e-3 1e-2 1e-1 1]);

set(ax, ...
    'XScale','log', ...
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

lgd = legend(ax,[hCurve;hTolerance], ...
    [legendText;{'$\epsilon_{\mathrm{rel}}=10^{-3}$'}], ...
    'Location','northeast', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex', ...
    'Box','on');

lgd.LineWidth = 1;
lgd.ItemTokenSize = [12 12];

%% ========================================================================
% 5. EXPORT FIGURE AND DATA
% ========================================================================

exportPublicationFigure(fig,'Fig09_Bound_Convergence.pdf','Fig09_Bound_Convergence.png');

SNR_dB = [];
Number_of_Intervals = [];
Lower_Bound_OP = [];
Upper_Bound_OP = [];
Absolute_Gap = [];
Relative_Gap = [];

for caseIndex = 1:numberOfCases

    numberOfRows = numel(history{caseIndex}.numberOfIntervals);

    SNR_dB = [SNR_dB; ...
        repmat(snrDbCases(caseIndex),numberOfRows,1)]; %#ok<AGROW>

    Number_of_Intervals = [Number_of_Intervals; ...
        history{caseIndex}.numberOfIntervals(:)]; %#ok<AGROW>

    Lower_Bound_OP = [Lower_Bound_OP; ...
        history{caseIndex}.lowerBound(:)]; %#ok<AGROW>

    Upper_Bound_OP = [Upper_Bound_OP; ...
        history{caseIndex}.upperBound(:)]; %#ok<AGROW>

    Absolute_Gap = [Absolute_Gap; ...
        history{caseIndex}.absoluteGap(:)]; %#ok<AGROW>

    Relative_Gap = [Relative_Gap; ...
        history{caseIndex}.relativeGap(:)]; %#ok<AGROW>
end

resultsTable = table( ...
    SNR_dB,Number_of_Intervals, ...
    Lower_Bound_OP,Upper_Bound_OP, ...
    Absolute_Gap,Relative_Gap);

writetable(resultsTable,'Fig09_Bound_Convergence_Data.csv');

save('Fig09_Bound_Convergence_Data.mat', ...
    'resultsTable','history','snrDbCases','snrLinCases', ...
    'targetRelativeGap','displayTolerance');

fprintf('\nExported:\n');
fprintf('Fig09_Bound_Convergence.pdf\n');
fprintf('Fig09_Bound_Convergence.png\n');
fprintf('Fig09_Bound_Convergence_Data.csv\n');
fprintf('Fig09_Bound_Convergence_Data.mat\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function result = adaptiveBoundHistoryProbabilityMidpoint( ...
    a,b,gammaTh,Ns,Nd,dist1,dist2, ...
    targetRelativeGap,maximumIntervals)
% Records the complete convergence history of the finite-sum bounds.
%
% At every iteration, the interval having the largest local finite-partition gap
% is divided at the midpoint of its probability mass:
%
%   F_X(x_mid) = [F_X(x_left)+F_X(x_right)]/2.
%
% This produces two subintervals with equal first-hop probability masses.

    x0 = gammaTh/a;
    F0 = selectedCDF(x0,dist1,Ns);

    % Initial partition [x0,infinity].
    xNodes = [x0,Inf];
    FNodes = [F0,1];

    % H(x) = F_Y(T(x)).
    % H(x0+) = 1 and H(infinity) = 0.
    HNodes = [1,0];

    currentIntervals = 1;

    intervalHistory    = [];
    lowerHistory       = [];
    upperHistory       = [];
    absoluteGapHistory = [];
    relativeGapHistory = [];

    while true

        probabilityMass = diff(FNodes);

        lowerContributions = ...
            HNodes(2:end).*probabilityMass;

        upperContributions = ...
            HNodes(1:end-1).*probabilityMass;

        P_LB = F0+sum(lowerContributions);
        P_UB = F0+sum(upperContributions);

        absoluteGap = P_UB-P_LB;
        midpointEstimate = 0.5*(P_UB+P_LB);

        relativeGap = ...
            absoluteGap/max(midpointEstimate,realmin);

        intervalHistory(end+1,1) = currentIntervals; %#ok<AGROW>
        lowerHistory(end+1,1) = P_LB; %#ok<AGROW>
        upperHistory(end+1,1) = P_UB; %#ok<AGROW>
        absoluteGapHistory(end+1,1) = absoluteGap; %#ok<AGROW>
        relativeGapHistory(end+1,1) = relativeGap; %#ok<AGROW>

        if relativeGap <= targetRelativeGap
            break;
        end

        if currentIntervals >= maximumIntervals
            warning(['Maximum number of intervals reached before the ' ...
                     'target relative gap was achieved.']);
            break;
        end

        % Finite-partition contribution of each interval to the total bound gap.
        localGap = ...
            (HNodes(1:end-1)-HNodes(2:end)).* ...
            probabilityMass;

        % Refine the interval with the largest local contribution.
        [~,intervalIndex] = max(localGap);

        xLeft  = xNodes(intervalIndex);
        xRight = xNodes(intervalIndex+1);

        FLeft  = FNodes(intervalIndex);
        FRight = FNodes(intervalIndex+1);

        % Midpoint of the interval probability mass.
        targetProbability = 0.5*(FLeft+FRight);

        xMid = inverseSelectedCDF( ...
            targetProbability,dist1,Ns,xLeft,xRight);

        FMid = selectedCDF(xMid,dist1,Ns);

        thresholdMid = exactConditionalThreshold( ...
            xMid,a,b,gammaTh);

        HMid = selectedCDF(thresholdMid,dist2,Nd);

        % Insert the new partition point.
        xNodes = [ ...
            xNodes(1:intervalIndex), ...
            xMid, ...
            xNodes(intervalIndex+1:end)];

        FNodes = [ ...
            FNodes(1:intervalIndex), ...
            FMid, ...
            FNodes(intervalIndex+1:end)];

        HNodes = [ ...
            HNodes(1:intervalIndex), ...
            HMid, ...
            HNodes(intervalIndex+1:end)];

        currentIntervals = currentIntervals+1;
    end

    result.numberOfIntervals = intervalHistory;
    result.lowerBound = lowerHistory;
    result.upperBound = upperHistory;
    result.absoluteGap = absoluteGapHistory;
    result.relativeGap = relativeGapHistory;
    result.finalXNodes = xNodes;
end

function xValue = inverseSelectedCDF( ...
    targetProbability,dist,numberOfPositions,xLower,xUpper)
% Solves F_X(xValue) = targetProbability using bisection.

    if targetProbability <= 0
        xValue = 0;
        return;
    end

    if targetProbability >= 1
        xValue = Inf;
        return;
    end

    % Construct a finite upper bracket when the interval ends at infinity.
    if isinf(xUpper)

        xUpper = max(1,2*xLower+1);

        while selectedCDF( ...
                xUpper,dist,numberOfPositions) < targetProbability

            xUpper = 2*xUpper+1;

            if xUpper > 1e12
                error('Unable to bracket the selected-CDF inverse.');
            end
        end
    end

    % Ensure the lower endpoint is finite.
    if ~isfinite(xLower)
        error('The lower inverse-CDF bracket must be finite.');
    end

    % Bisection inversion.
    for iteration = 1:80

        xMiddle = 0.5*(xLower+xUpper);

        FMiddle = selectedCDF( ...
            xMiddle,dist,numberOfPositions);

        if FMiddle < targetProbability
            xLower = xMiddle;
        else
            xUpper = xMiddle;
        end
    end

    xValue = 0.5*(xLower+xUpper);
end

function threshold = exactConditionalThreshold(x,a,b,gammaTh)

    denominator = b.*x.*(a.*x-gammaTh);

    threshold = inf(size(x));

    valid = denominator > 0;

    threshold(valid) = ...
        gammaTh.*(a.*x(valid)+1) ./ ...
        denominator(valid);
end

function dist = buildAlphaEtaMuDistribution( ...
    alpha,eta,mu,Omega,format,tolerance,maximumTerms)

    validateattributes(alpha,{'numeric'},{'scalar','positive'});
    validateattributes(mu,{'numeric'},{'scalar','positive'});
    validateattributes(Omega,{'numeric'},{'scalar','positive'});

    [h,H] = etaMuAuxiliaryParameters(eta,format);

    if abs(H) < 1e-15

        k = 0;
        weights = 1;

    else

        k = 0:maximumTerms;

        logWeights = ...
            gammaln(mu+k) ...
            -gammaln(mu) ...
            -gammaln(k+1) ...
            -mu*log(h) ...
            +k.*log((H/h)^2);

        weights = exp(logWeights);
        cumulativeWeights = cumsum(weights);

        finalIndex = find( ...
            cumulativeWeights >= 1-tolerance, ...
            1,'first');

        if isempty(finalIndex)

            warning(['The alpha-eta-mu series did not reach the ' ...
                     'requested tolerance. Increase maximumTerms.']);

            finalIndex = numel(k);
        end

        k = k(1:finalIndex);
        weights = weights(1:finalIndex);

        weights = weights/sum(weights);
    end

    dist.alpha = alpha;
    dist.eta = eta;
    dist.mu = mu;
    dist.Omega = Omega;
    dist.format = format;
    dist.h = h;
    dist.H = H;
    dist.k = k;
    dist.shape = 2*(mu+k);
    dist.weights = weights;
end

function [h,H] = etaMuAuxiliaryParameters(eta,format)

    switch format

        case 1
            if eta <= 0
                error('Format I requires eta > 0.');
            end

            h = (1+eta)^2/(4*eta);
            H = (1-eta^2)/(4*eta);

        case 2
            if eta <= -1 || eta >= 1
                error('Format II requires -1 < eta < 1.');
            end

            h = 1/(1-eta^2);
            H = eta/(1-eta^2);

        otherwise
            error('format must be either 1 or 2.');
    end
end

function F = selectedCDF(z,dist,numberOfPositions)

    F = parentCDF(z,dist).^numberOfPositions;
end

function F = parentCDF(z,dist)

    F = zeros(size(z));

    positiveFinite   = (z > 0) & isfinite(z);
    positiveInfinity = isinf(z) & (z > 0);

    if any(positiveFinite(:))

        localZ = z(positiveFinite);

        argument = ...
            2*dist.mu*dist.h .* ...
            (localZ./dist.Omega).^(dist.alpha/2);

        localCDF = zeros(size(localZ));

        for index = 1:numel(dist.weights)

            localCDF = localCDF + ...
                dist.weights(index).* ...
                gammainc( ...
                argument, ...
                dist.shape(index), ...
                'lower');
        end

        F(positiveFinite) = localCDF;
    end

    F(positiveInfinity) = 1;

    F = min(max(real(F),0),1);
end

function exportPublicationFigure(fig,pdfFile,pngFile)

    try
        exportgraphics( ...
            fig,pdfFile, ...
            'ContentType','vector', ...
            'BackgroundColor','white');

        exportgraphics( ...
            fig,pngFile, ...
            'Resolution',600, ...
            'BackgroundColor','white');

    catch
        set(fig,'PaperPositionMode','auto');
        print(fig,pdfFile,'-dpdf','-painters');
        print(fig,pngFile,'-dpng','-r600');
    end
end
