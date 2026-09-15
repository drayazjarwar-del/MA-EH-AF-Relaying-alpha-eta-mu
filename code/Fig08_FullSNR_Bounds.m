
% Publication color palette (same as Figures 2-5)
curveColors = [
    0.0000 0.4470 0.7410;   % Blue
    0.8500 0.3250 0.0980;   % Orange
    0.4660 0.6740 0.1880;   % Green
    0.4940 0.1840 0.5560    % Purple
];
%% Figure 8: Full-SNR finite-sum bounds with finite-partition error control for the PSR protocol
% MA-assisted energy-harvesting variable-gain AF relaying over
% alpha-eta-mu fading.
%
% This publication-ready script generates the Fig. 8 using the
% ORIGINAL numerical parameters of the paper:
%
%   alpha_1 = alpha_2 = 2
%   eta_1   = eta_2   = 0.5
%   mu_1    = mu_2    = 2
%   Ns = Nd = 3
%   zeta = 0.8
%   lambda = 0.6
%   R = 2 bit/s/Hz
%
% For PSR:
%
%   a = (1-lambda)*gammaBar,
%   b = zeta*lambda*gammaBar,
%   gamma_th = 2^(2R)-1.
%
% The figure contains:
%   - Finite-sum lower bound: dotted black line
%   - Exact analytical result from Eq. (24): solid black line
%   - Finite-sum upper bound: dashed black line
%   - Monte Carlo simulation: open black circles
%   - Inset: normalized lower/exact and upper/exact ratios
%
% Exported files:
%   Fig08_FullSNR_Bounds.pdf
%   Fig08_FullSNR_Bounds.png
%   Fig08_Bounds_Data.csv
%   Fig08_Bounds_Data.mat
%
% Public release filename:
%   Fig08_FullSNR_Bounds.m

clear;
clc;
close all;

rng(20260711,'twister');

%% ========================================================================
% 1. PAPER PARAMETERS
% ========================================================================

% First hop: S -> R
alpha1  = 2;
eta1    = 0.5;
mu1     = 2;
Omega1  = 1;
format1 = 1;                  % 1 = Format I, 2 = Format II
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
R      = 2;                  % bit/s/Hz

% Average source transmit SNR
snrDb  = 0:2:30;
snrLin = 10.^(snrDb/10);

% Physical PSR threshold used in the paper
gammaTh = 2^(2*R)-1;

% Protocol constants
A = 1-lambda;
B = zeta*lambda;

% Monte Carlo settings
Nmc = 1e6;

% Alpha-eta-mu series settings
seriesTolerance = 1e-14;
maximumTerms    = 500;

% Exact one-dimensional outage-integration settings
exactRelTol = 1e-9;
exactAbsTol = 1e-20;

% Finite-partition bound settings
% A relative tolerance of 10^(-3) gives approximately 0.1% separation.
boundRelTol      = 1e-3;
boundAbsTol      = 1e-14;
maximumIntervals = 5000;

%% ========================================================================
% 2. BUILD PARENT DISTRIBUTIONS
% ========================================================================

dist1 = buildAlphaEtaMuDistribution( ...
    alpha1,eta1,mu1,Omega1,format1, ...
    seriesTolerance,maximumTerms);

dist2 = buildAlphaEtaMuDistribution( ...
    alpha2,eta2,mu2,Omega2,format2, ...
    seriesTolerance,maximumTerms);

%% ========================================================================
% 3. GENERATE MONTE CARLO CHANNEL REALIZATIONS
% ========================================================================

fprintf('Generating %g Monte Carlo channel realizations...\n',Nmc);

Xmc = sampleSelectedAlphaEtaMu(Nmc,Ns,dist1);
Ymc = sampleSelectedAlphaEtaMu(Nmc,Nd,dist2);

fprintf('Channel generation completed.\n\n');

%% ========================================================================
% 4. EXACT, SIMULATION, AND FINITE-SUM BOUNDS
% ========================================================================

numberOfSnrPoints = numel(snrDb);

Pexact = zeros(1,numberOfSnrPoints);
Psim   = zeros(1,numberOfSnrPoints);
Plower = zeros(1,numberOfSnrPoints);
Pupper = zeros(1,numberOfSnrPoints);

outageCount    = zeros(1,numberOfSnrPoints);
boundGap       = zeros(1,numberOfSnrPoints);
relativeGap    = zeros(1,numberOfSnrPoints);
numberIntervals = zeros(1,numberOfSnrPoints);

for snrIndex = 1:numberOfSnrPoints

    gammaBar = snrLin(snrIndex);

    a = A*gammaBar;
    b = B*gammaBar;

    % Exact analytical outage probability from the one-dimensional representation
    Pexact(snrIndex) = exactOutageEq24( ...
        a,b,gammaTh,Ns,Nd,dist1,dist2, ...
        exactRelTol,exactAbsTol);

    % Direct Monte Carlo simulation
    gamma1 = a.*Xmc;
    gamma2 = b.*Xmc.*Ymc;

    gammaAF = (gamma1.*gamma2) ./ ...
              (gamma1+gamma2+1);

    outageCount(snrIndex) = sum(gammaAF < gammaTh);
    Psim(snrIndex) = outageCount(snrIndex)/Nmc;

    % Adaptive finite-sum bounds with finite-partition error control
    [Plower(snrIndex),Pupper(snrIndex),boundInformation] = ...
        adaptiveFullSnrBounds( ...
        a,b,gammaTh,Ns,Nd,dist1,dist2, ...
        boundRelTol,boundAbsTol,maximumIntervals);

    boundGap(snrIndex) = Pupper(snrIndex)-Plower(snrIndex);

    relativeGap(snrIndex) = ...
        boundGap(snrIndex) / ...
        max(0.5*(Pupper(snrIndex)+Plower(snrIndex)),realmin);

    numberIntervals(snrIndex) = ...
        boundInformation.numberOfIntervals;

    % Numerical ordering check
    orderingTolerance = max( ...
        10*exactAbsTol, ...
        1e-7*max(Pexact(snrIndex),1e-12));

    if Plower(snrIndex) > Pexact(snrIndex)+orderingTolerance
        warning(['Lower bound exceeds the exact result at %.1f dB. ' ...
                 'Tighten the numerical tolerances.'],snrDb(snrIndex));
    end

    if Pexact(snrIndex) > Pupper(snrIndex)+orderingTolerance
        warning(['Exact result exceeds the upper bound at %.1f dB. ' ...
                 'Tighten the numerical tolerances.'],snrDb(snrIndex));
    end

    fprintf(['SNR = %2.0f dB | Exact = %.6e | Sim = %.6e | ' ...
             'LB = %.6e | UB = %.6e | Rel. gap = %.3e | L = %d\n'], ...
             snrDb(snrIndex), ...
             Pexact(snrIndex), ...
             Psim(snrIndex), ...
             Plower(snrIndex), ...
             Pupper(snrIndex), ...
             relativeGap(snrIndex), ...
             numberIntervals(snrIndex));
end

%% ========================================================================
% 5. PUBLICATION-READY FIGURE
% ========================================================================

% Times font and compact single-column dimensions
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 3.5 3.10], ...
    'PaperPositionMode','auto');

mainAxes = axes(fig);
hold(mainAxes,'on');
box(mainAxes,'on');
grid(mainAxes,'on');

% Plot bounds first, then exact result, so the solid exact curve remains
% visible when all three curves are very close.
semilogy(mainAxes, ...
    snrDb,Plower,':', ...
    'Color',[0 0.4470 0.7410], ...
    'LineWidth',1.2, ...
    'DisplayName','Lower bound');

semilogy(mainAxes, ...
    snrDb,Pupper,'--', ...
    'Color',[0.8500 0.3250 0.0980], ...
    'LineWidth',1.2, ...
    'DisplayName','Upper bound');

semilogy(mainAxes, ...
    snrDb,Pexact,'-', ...
    'Color',[0.4660 0.6740 0.1880], ...
    'LineWidth',1.2, ...
    'DisplayName','Exact');

validSimulation = Psim > 0;

semilogy(mainAxes, ...
    snrDb(validSimulation), ...
    Psim(validSimulation),'o', ...
    'Color','k', ...
    'LineStyle','none', ...
    'MarkerSize',6, ...
    'LineWidth',0.8, ...
    'MarkerFaceColor','none', ...
    'DisplayName','Simulation');

xlabel(mainAxes,'Average source transmit SNR (dB)', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');
ylabel(mainAxes,'Outage probability', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

xlim(mainAxes,[0 30]);
xticks(mainAxes,0:5:30);

% Keep the same outage range used in the existing numerical figures.
ylim(mainAxes,[1e-6 1]);
yticks(mainAxes,10.^(-6:1:0));

set(mainAxes, ...
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

legend(mainAxes, ...
    'Location','southwest', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Box','on');

% ---------------- Normalized-bound inset ----------------
ratioLower = Plower./Pexact;
ratioUpper = Pupper./Pexact;

insetAxes = axes( ...
    fig, ...
    'Position',[0.3 0.25 0.32 0.26]);

hold(insetAxes,'on');
box(insetAxes,'on');
grid(insetAxes,'on');

plot(insetAxes, ...
    snrDb,ratioLower,':', ...
    'Color',[0 0.4470 0.7410], ...
    'LineWidth',1.0, ...
    'DisplayName','$P_{\rm out}^{\rm LB}/P_{\rm out}$');

plot(insetAxes, ...
    snrDb,ratioUpper,'--', ...
    'Color',[0.8500 0.3250 0.0980], ...
    'LineWidth',1.0, ...
    'DisplayName','$P_{\rm out}^{\rm UB}/P_{\rm out}$');

plot(insetAxes,[snrDb(1) snrDb(end)],[1 1],'-', ...
    'Color',[0.4660 0.6740 0.1880], ...
    'LineWidth',0.6, ...
    'HandleVisibility','off');

xlim(insetAxes,[10 30]);
xticks(insetAxes,10:10:30);

minimumRatio = min([ratioLower ratioUpper]);
maximumRatio = max([ratioLower ratioUpper]);

ratioRange = maximumRatio-minimumRatio;
ratioPadding = max(2e-4,0.20*ratioRange);

insetLower = min(minimumRatio-ratioPadding,1-ratioPadding);
insetUpper = max(maximumRatio+ratioPadding,1+ratioPadding);

ylim(insetAxes,[insetLower insetUpper]);

set(insetAxes, ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'LineWidth',0.65, ...
    'TickDir','in');

xlabel(insetAxes,'SNR (dB)', ...
    'FontName','Times New Roman', ...
    'FontSize',8);
ylabel(insetAxes,'Normalised bound', ...
    'FontName','Times New Roman', ...
    'FontSize',8);


%% ========================================================================
% 6. EXPORT FIGURE AND NUMERICAL DATA
% ========================================================================

exportPublicationFigure(fig,'Fig08_FullSNR_Bounds.pdf','Fig08_FullSNR_Bounds.png');

SNR_dB          = snrDb(:);
Exact_OP        = Pexact(:);
Simulation_OP   = Psim(:);
Lower_Bound_OP  = Plower(:);
Upper_Bound_OP  = Pupper(:);
Relative_Gap    = relativeGap(:);
Intervals       = numberIntervals(:);
Outage_Count    = outageCount(:);

resultsTable = table( ...
    SNR_dB,Exact_OP,Simulation_OP, ...
    Lower_Bound_OP,Upper_Bound_OP, ...
    Relative_Gap,Intervals,Outage_Count);

writetable(resultsTable,'Fig08_Bounds_Data.csv');

save('Fig08_Bounds_Data.mat', ...
    'resultsTable','snrDb','snrLin', ...
    'Pexact','Psim','Plower','Pupper', ...
    'relativeGap','numberIntervals','outageCount');

fprintf('\nExported:\n');
fprintf('Fig08_FullSNR_Bounds.pdf\n');
fprintf('Fig08_FullSNR_Bounds.png\n');
fprintf('Fig08_Bounds_Data.csv\n');
fprintf('Fig08_Bounds_Data.mat\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [P_LB,P_UB,information] = adaptiveFullSnrBounds( ...
    a,b,gammaTh,Ns,Nd,dist1,dist2, ...
    relativeTolerance,absoluteTolerance,maximumIntervals)
% Adaptive finite-partition-controlled bounds derived directly from the exact outage representation.

    x0 = gammaTh/a;
    F0 = selectedCDF(x0,dist1,Ns);

    % Initial transformed partition:
    % t = 0 -> x = x0, t = 1 -> x = infinity.
    tNodes = [0,1];
    xNodes = [x0,Inf];

    FNodes = [F0,1];

    % H(x) = F_Y(T(x)); H(x0+) = 1 and H(infinity) = 0.
    HNodes = [1,0];

    numberOfIntervals = 1;
    converged = false;

    while true

        probabilityMass = diff(FNodes);

        lowerContributions = ...
            HNodes(2:end).*probabilityMass;

        upperContributions = ...
            HNodes(1:end-1).*probabilityMass;

        P_LB = F0+sum(lowerContributions);
        P_UB = F0+sum(upperContributions);

        localGap = ...
            (HNodes(1:end-1)-HNodes(2:end)).* ...
            probabilityMass;

        totalGap = P_UB-P_LB;
        midpointEstimate = 0.5*(P_LB+P_UB);

        stoppingTolerance = max( ...
            absoluteTolerance, ...
            relativeTolerance*max(midpointEstimate,realmin));

        if totalGap <= stoppingTolerance
            converged = true;
            break;
        end

        if numberOfIntervals >= maximumIntervals
            warning(['Maximum number of partition intervals reached ' ...
                     'before satisfying the requested bound tolerance.']);
            break;
        end

        % Bisect the interval contributing the largest finite-partition gap.
        [~,intervalIndex] = max(localGap);

        tLeft  = tNodes(intervalIndex);
        tRight = tNodes(intervalIndex+1);
        tMid   = 0.5*(tLeft+tRight);

        xMid = x0+tMid/(1-tMid);

        FMid = selectedCDF(xMid,dist1,Ns);

        thresholdMid = exactConditionalThreshold( ...
            xMid,a,b,gammaTh);

        HMid = selectedCDF(thresholdMid,dist2,Nd);

        tNodes = [ ...
            tNodes(1:intervalIndex), ...
            tMid, ...
            tNodes(intervalIndex+1:end)];

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

        numberOfIntervals = numberOfIntervals+1;
    end

    P_LB = min(max(real(P_LB),0),1);
    P_UB = min(max(real(P_UB),0),1);

    information.numberOfIntervals = numberOfIntervals;
    information.finitePartitionGap = P_UB-P_LB;
    information.converged = converged;
    information.tNodes = tNodes;
    information.xNodes = xNodes;
end

function threshold = exactConditionalThreshold(x,a,b,gammaTh)

    denominator = b.*x.*(a.*x-gammaTh);

    threshold = inf(size(x));

    valid = denominator > 0;

    threshold(valid) = ...
        gammaTh.*(a.*x(valid)+1) ./ ...
        denominator(valid);
end

function Pout = exactOutageEq24( ...
    a,b,gammaTh,Ns,Nd,dist1,dist2, ...
    relativeTolerance,absoluteTolerance)

    x0 = gammaTh/a;

    firstTerm = selectedCDF(x0,dist1,Ns);

    integrand = @(t) exactEq24TransformedIntegrand( ...
        t,x0,a,b,gammaTh,Ns,Nd,dist1,dist2);

    integralPartA = integral( ...
        integrand,0,0.9, ...
        'RelTol',relativeTolerance, ...
        'AbsTol',absoluteTolerance);

    integralPartB = integral( ...
        integrand,0.9,1, ...
        'RelTol',relativeTolerance, ...
        'AbsTol',absoluteTolerance);

    Pout = firstTerm+integralPartA+integralPartB;

    Pout = min(max(real(Pout),0),1);
end

function value = exactEq24TransformedIntegrand( ...
    t,x0,a,b,gammaTh,Ns,Nd,dist1,dist2)

    value = zeros(size(t));

    valid = t < 1;

    if ~any(valid(:))
        return;
    end

    localT = t(valid);

    x = x0+localT./(1-localT);
    dxDt = 1./(1-localT).^2;

    yThreshold = exactConditionalThreshold( ...
        x,a,b,gammaTh);

    FY = selectedCDF(yThreshold,dist2,Nd);
    fX = selectedPDF(x,dist1,Ns);

    localValue = FY.*fX.*dxDt;
    localValue(~isfinite(localValue)) = 0;

    value(valid) = localValue;
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
    dist.cumulativeWeights = cumsum(weights);
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

function selectedGain = sampleSelectedAlphaEtaMu( ...
    numberOfSamples,numberOfPositions,dist)

    selectedGain = zeros(numberOfSamples,1);

    for positionIndex = 1:numberOfPositions

        candidateGain = sampleParentAlphaEtaMu( ...
            numberOfSamples,dist);

        selectedGain = max(selectedGain,candidateGain);
    end
end

function Z = sampleParentAlphaEtaMu(numberOfSamples,dist)

    uniformRandom = rand(numberOfSamples,1);

    mixtureIndex = ones(numberOfSamples,1);

    for index = 1:numel(dist.weights)-1

        mixtureIndex( ...
            uniformRandom > dist.cumulativeWeights(index)) = ...
            index+1;
    end

    shape = dist.shape(mixtureIndex);

    U = randg(shape(:));

    Z = dist.Omega .* ...
        (U./(2*dist.mu*dist.h)).^(2/dist.alpha);
end

function F = selectedCDF(z,dist,numberOfPositions)

    F = parentCDF(z,dist).^numberOfPositions;
end

function f = selectedPDF(z,dist,numberOfPositions)

    Fparent = parentCDF(z,dist);
    fparent = parentPDF(z,dist);

    f = numberOfPositions.*fparent.* ...
        Fparent.^(numberOfPositions-1);
end

function F = parentCDF(z,dist)

    F = zeros(size(z));

    positiveFinite = (z > 0) & isfinite(z);
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

function f = parentPDF(z,dist)

    f = zeros(size(z));

    positiveFinite = (z > 0) & isfinite(z);

    if any(positiveFinite(:))

        localZ = z(positiveFinite);

        argument = ...
            2*dist.mu*dist.h .* ...
            (localZ./dist.Omega).^(dist.alpha/2);

        derivative = ...
            (dist.alpha*dist.mu*dist.h / ...
            dist.Omega^(dist.alpha/2)) .* ...
            localZ.^(dist.alpha/2-1);

        mixturePDF = zeros(size(localZ));

        for index = 1:numel(dist.weights)

            shape = dist.shape(index);

            gammaDensity = exp( ...
                (shape-1).*log(argument) ...
                -argument ...
                -gammaln(shape));

            mixturePDF = mixturePDF + ...
                dist.weights(index).*gammaDensity;
        end

        f(positiveFinite) = ...
            mixturePDF.*derivative;
    end

    f(~isfinite(f)) = 0;
    f = max(real(f),0);
end

function exportPublicationFigure(fig,pdfFile,pngFile)

    try
        print(fig,pdfFile,'-dpdf','-painters');

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