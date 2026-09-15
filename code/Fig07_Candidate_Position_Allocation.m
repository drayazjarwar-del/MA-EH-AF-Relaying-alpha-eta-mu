%% FIGURE 7: Required SNR under optimal MA-position allocation
% Movable-antenna-assisted energy-harvesting AF relaying over
% alpha-eta-mu fading channels.
%
% DESIGN QUESTION
% ---------------
% Given a fixed total candidate-position budget
%
%       N_T = N_s + N_d = 8,
%
% how should the candidate positions be divided between the source and
% destination to minimize the average source transmit SNR required to
% satisfy
%
%       P_out <= 10^(-3)?
%
% The calculation is repeated for:
%
%       lambda = 0.2, 0.6, and 0.8
%
% under the PSR protocol.
%
% ANALYTICAL RESULT
% -----------------
% For every (N_s,N_d) pair, the exact one-dimensional outage representation is solved with respect to average source transmit SNR.
%
% MONTE CARLO RESULT
% ------------------
% For each realization, the exact source SNR required to satisfy the AF
% threshold is calculated from a quadratic equation. A streaming histogram
% then estimates its (1-P_target) quantile without storing 10^7 samples.
%
% OUTPUTS
% -------
%   Fig07_Candidate_Position_Allocation.pdf       vector PDF for LaTeX
%   Fig07_Candidate_Position_Allocation.png       600-dpi backup
%   Fig07_Candidate_Position_Allocation.eps
%   Fig07_Candidate_Position_Allocation_results.csv
%
% MATLAB requirement:
%   Local functions in scripts require MATLAB R2016b or newer.
%
% No Statistics and Machine Learning Toolbox is required.

clear; clc; close all;
rng(2026,'twister');

set(groot,'defaultFigureVisible','on');

%% ------------------------ Visible progress window -----------------------
statusFig = figure('Color','w', ...
                   'Visible','on', ...
                   'Name','Figure 7 calculation status', ...
                   'NumberTitle','off', ...
                   'Position',[300 300 580 200]);

statusAx = axes(statusFig);
axis(statusAx,'off');

statusText = text(statusAx,0.5,0.61, ...
    'Figure 7 calculation is running...', ...
    'HorizontalAlignment','center', ...
    'FontSize',10);

statusSubText = text(statusAx,0.5,0.37, ...
    'The exact outage equation is being solved for the required SNR.', ...
    'HorizontalAlignment','center', ...
    'FontSize',10);

drawnow;

%% ---------------------------- Configuration -----------------------------
cfg.alpha1 = 2;
cfg.eta1   = 0.5;
cfg.mu1    = 2;
cfg.Omega1 = 1;
cfg.format1 = 1;                    % Format I

cfg.alpha2 = 2;
cfg.eta2   = 0.5;
cfg.mu2    = 2;
cfg.Omega2 = 1;
cfg.format2 = 1;                    % Format I

cfg.zeta = 0.8;
cfg.R    = 2;                       % target rate, bit/s/Hz

cfg.lambdaValues = [0.2,0.6,0.8];

cfg.NT = 8;                         % fixed total candidate-position budget
cfg.NsValues = 1:(cfg.NT-1);
cfg.NdValues = cfg.NT-cfg.NsValues;

cfg.targetOutage = 1e-3;

% Exact analytical root-search settings
cfg.snrLower_dB = -10;
cfg.snrUpper_dB = 50;
cfg.snrTolerance_dB = 1e-3;
cfg.maxBisectionIterations = 30;

% Parent-distribution and numerical-integral controls
cfg.seriesTol = 1e-13;
cfg.Kmax      = 500;
cfg.RelTol    = 1e-7;
cfg.AbsTol    = 1e-11;

% Monte Carlo controls
cfg.Nmc = 1e8;                      % final journal result
% For a quick test, temporarily use:
% cfg.Nmc = 1e5;

cfg.chunkSize = 1e5;

% Streaming histogram used to estimate the 99.9th percentile of the
% per-realization required SNR.
cfg.histMin_dB  = -20;
cfg.histMax_dB  = 100;
cfg.histStep_dB = 0.01;
cfg.histEdges_dB = cfg.histMin_dB:cfg.histStep_dB:cfg.histMax_dB;

% Plot styles
lineStyles = {'-','--','-.'};
lambdaLabels = { ...
    '$\lambda=0.2$', ...
    '$\lambda=0.6$', ...
    '$\lambda=0.8$'};

nLambda = numel(cfg.lambdaValues);
nAlloc  = numel(cfg.NsValues);
nBins   = numel(cfg.histEdges_dB)-1;

requiredSNR_Ana_dB = zeros(nAlloc,nLambda);
requiredSNR_MC_dB  = zeros(nAlloc,nLambda);

%% ---------------- Build parent fading distributions ---------------------
fprintf('Building alpha-eta-mu parent distributions...\n');
statusText.String = 'Building parent fading distributions...';
drawnow;

mix1 = buildAEMMixture( ...
    cfg.alpha1,cfg.eta1,cfg.mu1,cfg.Omega1,cfg.format1, ...
    cfg.seriesTol,cfg.Kmax);

mix2 = buildAEMMixture( ...
    cfg.alpha2,cfg.eta2,cfg.mu2,cfg.Omega2,cfg.format2, ...
    cfg.seriesTol,cfg.Kmax);

pdfMass1 = integral(@(z) aemPDF(z,mix1),0,Inf, ...
                    'RelTol',1e-8,'AbsTol',1e-11);

pdfMass2 = integral(@(z) aemPDF(z,mix2),0,Inf, ...
                    'RelTol',1e-8,'AbsTol',1e-11);

fprintf('First-hop parent PDF mass : %.12f\n',pdfMass1);
fprintf('Second-hop parent PDF mass: %.12f\n\n',pdfMass2);

if abs(pdfMass1-1)>1e-6 || abs(pdfMass2-1)>1e-6
    warning(['A parent PDF does not integrate sufficiently close to one. ', ...
             'Increase cfg.Kmax or inspect the fading parameters.']);
end

%% ------------------- Exact analytical required SNR ----------------------
fprintf('Solving the exact outage equation for required SNR...\n');

for ell = 1:nLambda
    lambda = cfg.lambdaValues(ell);

    for k = 1:nAlloc
        Ns = cfg.NsValues(k);
        Nd = cfg.NdValues(k);

        statusText.String = sprintf( ...
            'Analytical calculation: \\lambda=%.1f, (N_s,N_d)=(%d,%d)', ...
            lambda,Ns,Nd);

        statusSubText.String = sprintf( ...
            'Allocation %d of %d; power-splitting case %d of %d', ...
            k,nAlloc,ell,nLambda);

        drawnow;

        requiredSNR_Ana_dB(k,ell) = requiredSNRExact( ...
            cfg.targetOutage,lambda,Ns,Nd,cfg,mix1,mix2);

        fprintf(['lambda=%.1f | (Ns,Nd)=(%d,%d) | ', ...
                 'Analytical required SNR=%.4f dB\n'], ...
                 lambda,Ns,Nd,requiredSNR_Ana_dB(k,ell));
    end
end

%% --------------- Monte Carlo required-SNR calculation ------------------
% For a fixed realization (X,Y) and a source transmit SNR s = Ps/N0:
%
%   gamma1 = A*s*X,             A = 1-lambda
%   gamma2 = B*s*X*Y,           B = zeta*lambda
%
% Solving gamma_AF = gamma_th for s gives the positive quadratic root:
%
% s_req = [D + sqrt(D^2 + 4*A*B*X^2*Y*gamma_th)]
%         / [2*A*B*X^2*Y],
%
% where D = gamma_th*X*(A+B*Y).
%
% P_out(s) = Pr(s_req > s). Therefore, the required SNR for target outage
% epsilon is the (1-epsilon) quantile of s_req.

fprintf('\nRunning Monte Carlo required-SNR calculation with %.0e realizations...\n', ...
        cfg.Nmc);

statusText.String = 'Running Monte Carlo required-SNR calculation...';
statusSubText.String = 'Generating movable-antenna channel realizations.';
drawnow;

histCounts = zeros(nBins,nAlloc,nLambda);
underflowCounts = zeros(nAlloc,nLambda);
overflowCounts  = zeros(nAlloc,nLambda);

gammaTh = 2^(2*cfg.R)-1;
processed = 0;
Nmax = cfg.NT-1;

while processed < cfg.Nmc
    thisChunk = min(cfg.chunkSize,cfg.Nmc-processed);

    % Generate enough independent candidate gains for all possible
    % allocations. Cumulative maxima provide the selected gain for each N.
    Xcand = zeros(thisChunk,Nmax);
    Ycand = zeros(thisChunk,Nmax);

    for n = 1:Nmax
        Xcand(:,n) = sampleAEM(mix1,thisChunk);
        Ycand(:,n) = sampleAEM(mix2,thisChunk);
    end

    Xcum = cummax(Xcand,2);
    Ycum = cummax(Ycand,2);

    for ell = 1:nLambda
        lambda = cfg.lambdaValues(ell);
        A = 1-lambda;
        B = cfg.zeta*lambda;

        for k = 1:nAlloc
            Ns = cfg.NsValues(k);
            Nd = cfg.NdValues(k);

            X = Xcum(:,Ns);
            Y = Ycum(:,Nd);

            quadraticCoefficient = A*B.*X.^2.*Y;
            linearMagnitude = gammaTh.*X.*(A+B.*Y);

            sRequired = (linearMagnitude + ...
                sqrt(linearMagnitude.^2 + ...
                4.*quadraticCoefficient.*gammaTh)) ./ ...
                (2.*quadraticCoefficient);

            sRequired_dB = 10.*log10(sRequired);

            underflowCounts(k,ell) = underflowCounts(k,ell) + ...
                sum(sRequired_dB < cfg.histMin_dB);

            overflowCounts(k,ell) = overflowCounts(k,ell) + ...
                sum(sRequired_dB >= cfg.histMax_dB);

            histCounts(:,k,ell) = histCounts(:,k,ell) + ...
                histcounts(sRequired_dB,cfg.histEdges_dB).';
        end
    end

    processed = processed + thisChunk;

    if mod(processed,1e6)==0 || processed==cfg.Nmc
        fprintf('  Processed %.0f of %.0f realizations\n', ...
                processed,cfg.Nmc);

        if isgraphics(statusText)
            statusText.String = sprintf( ...
                'Monte Carlo calculation: %.0f of %.0f realizations', ...
                processed,cfg.Nmc);

            statusSubText.String = ...
                'Estimating the 99.9th percentile of required source SNR.';
            drawnow;
        end
    end
end

% Convert each streaming histogram into the (1-targetOutage) quantile.
quantileProbability = 1-cfg.targetOutage;

for ell = 1:nLambda
    for k = 1:nAlloc
        requiredSNR_MC_dB(k,ell) = quantileFromHistogram( ...
            histCounts(:,k,ell), ...
            underflowCounts(k,ell), ...
            overflowCounts(k,ell), ...
            cfg.histEdges_dB, ...
            cfg.Nmc, ...
            quantileProbability);
    end
end

%% ---------------------- Report optimum allocations ----------------------
fprintf('\n============================================================\n');
fprintf('OPTIMAL CANDIDATE-POSITION ALLOCATIONS\n');
fprintf('Target outage probability: %.1e\n',cfg.targetOutage);
fprintf('Fixed total budget: NT=%d\n',cfg.NT);
fprintf('============================================================\n');

balancedIndex = find(cfg.NsValues==cfg.NT/2,1);
sourceHeavyIndex = find(cfg.NsValues==cfg.NT-1,1);

optimumIndex = zeros(nLambda,1);
optimumNs = zeros(nLambda,1);
optimumNd = zeros(nLambda,1);
optimumSNR = zeros(nLambda,1);
savingBalanced = zeros(nLambda,1);
savingSourceHeavy = zeros(nLambda,1);

for ell = 1:nLambda
    [optimumSNR(ell),optimumIndex(ell)] = ...
        min(requiredSNR_Ana_dB(:,ell));

    optimumNs(ell) = cfg.NsValues(optimumIndex(ell));
    optimumNd(ell) = cfg.NdValues(optimumIndex(ell));

    if ~isempty(balancedIndex)
        savingBalanced(ell) = ...
            requiredSNR_Ana_dB(balancedIndex,ell)-optimumSNR(ell);
    else
        savingBalanced(ell) = NaN;
    end

    if ~isempty(sourceHeavyIndex)
        savingSourceHeavy(ell) = ...
            requiredSNR_Ana_dB(sourceHeavyIndex,ell)-optimumSNR(ell);
    else
        savingSourceHeavy(ell) = NaN;
    end

    fprintf('\nlambda=%.1f\n',cfg.lambdaValues(ell));
    fprintf('  Analytical optimum: (Ns,Nd)=(%d,%d)\n', ...
            optimumNs(ell),optimumNd(ell));
    fprintf('  Minimum required SNR: %.4f dB\n',optimumSNR(ell));
    fprintf('  SNR saving over balanced allocation: %.4f dB\n', ...
            savingBalanced(ell));
    fprintf('  SNR saving over source-heavy allocation (%d,1): %.4f dB\n', ...
            cfg.NT-1,savingSourceHeavy(ell));

    fprintf('  Analytical versus Monte Carlo values:\n');

    for k = 1:nAlloc
        fprintf(['    (Ns,Nd)=(%d,%d): Analytical=%.4f dB, ', ...
                 'MC=%.4f dB\n'], ...
                 cfg.NsValues(k),cfg.NdValues(k), ...
                 requiredSNR_Ana_dB(k,ell), ...
                 requiredSNR_MC_dB(k,ell));
    end
end

%% ------------------------------ Plot -----------------------------------
if isgraphics(statusFig)
    close(statusFig);
end

% IEEE Transactions single-column formatting.
fig = figure('Color','w', ...
             'Units','inches', ...
             'Position',[1 1 3.5 3.10], ...
             'Visible','on');

layout = tiledlayout(fig,1,1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax = nexttile(layout);
hold(ax,'on');

hLine = gobjects(nLambda,1);
x = 1:nAlloc;

for ell = 1:nLambda
   
    
     % Monte Carlo required-SNR markers: common black open circles.
    plot(ax,x,requiredSNR_MC_dB(:,ell),'o', ...
        'LineStyle','none', ...
        'Color','k', ...
        'MarkerFaceColor','w', ...
        'MarkerSize',5, ...
        'LineWidth',0.6, ...
        'HandleVisibility','off');
    
    % Analytical curves: Figure 2/3 publication color style.
    curveColors = [0 0.4470 0.7410; ...
                   0.8500 0.3250 0.0980; ...
                   0.4660 0.6740 0.1880];

    hLine(ell) = plot(ax,x,requiredSNR_Ana_dB(:,ell), ...
        lineStyles{ell}, ...
        'Color',curveColors(ell,:), ...
        'LineWidth',1.2, ...
        'DisplayName',lambdaLabels{ell});

   
end

% Dummy marker for the Monte Carlo legend entry.
hSim = plot(ax,NaN,NaN,'ko', ...
    'LineStyle','none', ...
    'MarkerFaceColor','w', ...
    'MarkerSize',5, ...
    'LineWidth',0.6, ...
    'DisplayName','Simulation');

grid(ax,'on');
box(ax,'on');

ax.XMinorGrid = 'off';
ax.YMinorGrid = 'on';
ax.GridAlpha = 0.22;
ax.MinorGridAlpha = 0.12;
ax.FontName = 'Times New Roman';
ax.FontSize = 10;
ax.LineWidth = 0.8;
ax.TickDir = 'in';
ax.Layer = 'top';

allocationLabels = arrayfun( ...
    @(Ns,Nd) sprintf('(%d,%d)',Ns,Nd), ...
    cfg.NsValues,cfg.NdValues, ...
    'UniformOutput',false);

xticks(ax,x);
xticklabels(ax,allocationLabels);

xlabel(ax, ...
    'Candidate-position allocation $(N_s,N_d)$, $N_T=N_s+N_d=8$', ...
    'Interpreter','latex', ...
    'FontSize',10);

ylabel(ax, ...
    'Required average source transmit SNR (dB)', ...
    'Interpreter','latex', ...
    'FontSize',10);

% Vertical legend so that every entry appears on a separate line.
lgd = legend(ax,[hLine;hSim], ...
    [lambdaLabels,{'Simulation'}], ...
    'Location','northeast', ...
    'Orientation','vertical', ...
    'Interpreter','latex');

lgd.FontName = 'Times New Roman';
lgd.FontSize = 10;
lgd.Box = 'on';
lgd.LineWidth = 1;
lgd.ItemTokenSize = [12 12];

xlim(ax,[0.8 nAlloc+0.2]);

allValues = [requiredSNR_Ana_dB(:);requiredSNR_MC_dB(:)];
finiteValues = allValues(isfinite(allValues));

if ~isempty(finiteValues)
    yPadding = max(0.5,0.08*(max(finiteValues)-min(finiteValues)));
    ylim(ax,[min(finiteValues)-yPadding,max(finiteValues)+yPadding]);
end

drawnow;
figure(fig);
shg;

%% ---------------------------- Save outputs ------------------------------
print(fig,'Fig07_Candidate_Position_Allocation.pdf','-dpdf','-painters');

exportgraphics(fig,'Fig07_Candidate_Position_Allocation.png', ...
               'Resolution',600, ...
               'BackgroundColor','white');

print(fig,'Fig07_Candidate_Position_Allocation.eps','-depsc2','-r600');

T = table( ...
    cfg.NsValues(:),cfg.NdValues(:), ...
    requiredSNR_Ana_dB(:,1),requiredSNR_MC_dB(:,1), ...
    requiredSNR_Ana_dB(:,2),requiredSNR_MC_dB(:,2), ...
    requiredSNR_Ana_dB(:,3),requiredSNR_MC_dB(:,3), ...
    'VariableNames',{ ...
        'Ns','Nd', ...
        'lambda02_Analytical_dB','lambda02_MC_dB', ...
        'lambda06_Analytical_dB','lambda06_MC_dB', ...
        'lambda08_Analytical_dB','lambda08_MC_dB'});

writetable(T,'Fig07_Candidate_Position_Allocation_results.csv');

summaryTable = table( ...
    cfg.lambdaValues(:), ...
    optimumNs,optimumNd,optimumSNR, ...
    savingBalanced,savingSourceHeavy, ...
    'VariableNames',{ ...
        'lambda','Optimal_Ns','Optimal_Nd', ...
        'MinimumRequiredSNR_dB', ...
        'SavingVsBalanced_dB', ...
        'SavingVsSourceHeavy_dB'});

writetable(summaryTable,'Fig07_Optimal_Allocation_Summary.csv');

fprintf('\nFigure 4 and optimization results have been exported.\n');

%% =======================================================================
%                           LOCAL FUNCTIONS
% ========================================================================

function requiredSNR_dB = requiredSNRExact( ...
    targetOutage,lambda,Ns,Nd,cfg,mix1,mix2)
%REQUIREDSNREXACT Solve Pout(SNR)=targetOutage by bisection in dB.

    low_dB = cfg.snrLower_dB;
    high_dB = cfg.snrUpper_dB;

    pLow = outageAtSNRdB( ...
        low_dB,lambda,Ns,Nd,cfg,mix1,mix2);

    pHigh = outageAtSNRdB( ...
        high_dB,lambda,Ns,Nd,cfg,mix1,mix2);

    % Ensure the lower point is still in outage.
    lowerExpansionCount = 0;

    while pLow < targetOutage
        low_dB = low_dB-10;
        pLow = outageAtSNRdB( ...
            low_dB,lambda,Ns,Nd,cfg,mix1,mix2);

        lowerExpansionCount = lowerExpansionCount+1;

        if lowerExpansionCount>10
            error('Unable to bracket the target outage at low SNR.');
        end
    end

    % Ensure the upper point satisfies the outage requirement.
    upperExpansionCount = 0;

    while pHigh > targetOutage
        high_dB = high_dB+10;
        pHigh = outageAtSNRdB( ...
            high_dB,lambda,Ns,Nd,cfg,mix1,mix2);

        upperExpansionCount = upperExpansionCount+1;

        if upperExpansionCount>10
            error('Unable to bracket the target outage at high SNR.');
        end
    end

    iteration = 0;

    while (high_dB-low_dB)>cfg.snrTolerance_dB && ...
          iteration<cfg.maxBisectionIterations

        middle_dB = 0.5*(low_dB+high_dB);

        pMiddle = outageAtSNRdB( ...
            middle_dB,lambda,Ns,Nd,cfg,mix1,mix2);

        if pMiddle > targetOutage
            low_dB = middle_dB;
        else
            high_dB = middle_dB;
        end

        iteration = iteration+1;
    end

    requiredSNR_dB = 0.5*(low_dB+high_dB);
end

function Pout = outageAtSNRdB(snr_dB,lambda,Ns,Nd,cfg,mix1,mix2)
%OUTAGEATSNRDB Exact numerical outage under PSR at a specified SNR.

    snrLin = 10^(snr_dB/10);

    a = (1-lambda)*snrLin;
    b = cfg.zeta*lambda*snrLin;
    gammaTh = 2^(2*cfg.R)-1;

    Pout = outageExactNumerical( ...
        a,b,gammaTh,Ns,Nd,mix1,mix2,cfg.RelTol,cfg.AbsTol);
end

function q_dB = quantileFromHistogram( ...
    counts,underflowCount,overflowCount,edges_dB,totalSamples,qProbability)
%QUANTILEFROMHISTOGRAM Estimate a quantile from a streaming histogram.

    targetRank = qProbability*totalSamples;

    cumulativeBeforeBins = underflowCount;
    cumulativeBins = cumulativeBeforeBins+cumsum(counts);

    binIndex = find(cumulativeBins>=targetRank,1,'first');

    if isempty(binIndex)
        if overflowCount>0
            error(['The requested Monte Carlo quantile lies above ', ...
                   'cfg.histMax_dB. Increase cfg.histMax_dB.']);
        else
            error('Unable to determine the Monte Carlo quantile.');
        end
    end

    if binIndex==1
        countBefore = cumulativeBeforeBins;
    else
        countBefore = cumulativeBins(binIndex-1);
    end

    countInBin = counts(binIndex);

    if countInBin<=0
        q_dB = 0.5*(edges_dB(binIndex)+edges_dB(binIndex+1));
        return;
    end

    fractionWithinBin = ...
        (targetRank-countBefore)/countInBin;

    fractionWithinBin = min(max(fractionWithinBin,0),1);

    q_dB = edges_dB(binIndex) + ...
           fractionWithinBin* ...
           (edges_dB(binIndex+1)-edges_dB(binIndex));
end

function mix = buildAEMMixture(alpha,eta,mu,Omega,formatNo,tol,Kmax)
%BUILDAEMMIXTURE Gamma-mixture representation of alpha-eta-mu fading.

    validateattributes(alpha,{'numeric'},{'scalar','real','positive'});
    validateattributes(mu,{'numeric'},{'scalar','real','positive'});
    validateattributes(Omega,{'numeric'},{'scalar','real','positive'});

    [h,H] = etaAuxiliary(eta,formatNo);
    c = 2*mu*h;

    if abs(H)<1e-15
        k = 0;
        w = 1;
    else
        k = (0:Kmax).';

        logw = log(2)+0.5*log(pi) ...
             +2*mu*log(mu)+mu*log(h)-gammaln(mu) ...
             +2*k.*log(abs(mu*H)) ...
             -gammaln(k+1) ...
             -gammaln(mu+k+0.5) ...
             -2*(mu+k).*log(c) ...
             +gammaln(2*(mu+k));

        wAll = exp(logw);
        totalMass = sum(wAll);

        if abs(totalMass-1)>1e-8
            warning(['The series mass at Kmax=%d is %.12g rather than 1. ', ...
                     'Increase Kmax for these fading parameters.'], ...
                     Kmax,totalMass);
        end

        cumulativeMass = cumsum(wAll);
        last = find(cumulativeMass>=1-tol,1,'first');

        if isempty(last)
            last = numel(wAll);
        end

        k = k(1:last);
        w = wAll(1:last);
        w = w/sum(w);
    end

    mix.alpha = alpha;
    mix.eta = eta;
    mix.mu = mu;
    mix.Omega = Omega;
    mix.formatNo = formatNo;
    mix.h = h;
    mix.H = H;
    mix.c = c;
    mix.k = k;
    mix.shape = 2*(mu+k);
    mix.w = w;
end

function [h,H] = etaAuxiliary(eta,formatNo)
%ETAAUXILIARY Auxiliary h and H parameters.

    switch formatNo
        case 1
            if ~(eta>0)
                error('For Format I, eta must satisfy eta>0.');
            end

            h = (1+eta)^2/(4*eta);
            H = (1-eta^2)/(4*eta);

        case 2
            if ~(eta>-1 && eta<1)
                error('For Format II, eta must satisfy -1<eta<1.');
            end

            h = 1/(1-eta^2);
            H = eta/(1-eta^2);

        otherwise
            error('formatNo must be 1 or 2.');
    end
end

function F = aemCDF(z,mix)
%AEMCDF Parent alpha-eta-mu power-gain CDF.

    originalSize = size(z);
    z = z(:).';

    F = zeros(size(z));
    positive = z>0 & isfinite(z);

    if any(positive)
        u = mix.c.*(z(positive)./mix.Omega).^(mix.alpha/2);
        Fpositive = zeros(size(u));

        for k = 1:numel(mix.w)
            Fpositive = Fpositive + ...
                mix.w(k).*gammainc(u,mix.shape(k),'lower');
        end

        F(positive) = Fpositive;
    end

    F(isinf(z) & z>0) = 1;
    F(z<=0) = 0;
    F = min(max(F,0),1);
    F = reshape(F,originalSize);
end

function f = aemPDF(z,mix)
%AEMPDF Parent alpha-eta-mu power-gain PDF.

    originalSize = size(z);
    z = z(:).';

    f = zeros(size(z));
    positive = z>0 & isfinite(z);

    if any(positive)
        zPositive = z(positive);
        u = mix.c.*(zPositive./mix.Omega).^(mix.alpha/2);

        duDz = mix.c.*(mix.alpha/2)./mix.Omega ...
             .*(zPositive./mix.Omega).^(mix.alpha/2-1);

        fPositive = zeros(size(zPositive));

        for k = 1:numel(mix.w)
            shape = mix.shape(k);

            gammaPDF = exp( ...
                (shape-1).*log(u)-u-gammaln(shape));

            fPositive = fPositive + ...
                mix.w(k).*gammaPDF.*duDz;
        end

        f(positive) = fPositive;
    end

    f(z<=0) = 0;
    f = reshape(f,originalSize);
end

function z = sampleAEM(mix,N)
%SAMPLEAEM Generate N alpha-eta-mu parent-gain samples.

    uniformCategorical = rand(N,1);
    cumulativeWeights = cumsum(mix.w);
    cumulativeWeights(end) = 1;

    mixtureIndex = ones(N,1,'uint16');

    for k = 1:numel(cumulativeWeights)-1
        mixtureIndex(uniformCategorical>cumulativeWeights(k)) = k+1;
    end

    shapePerSample = mix.shape(double(mixtureIndex));
    gammaSample = randg(shapePerSample);

    z = mix.Omega.*(gammaSample./mix.c).^(2/mix.alpha);
end

function Pout = outageExactNumerical( ...
    a,b,gammaTh,Ns,Nd,mix1,mix2,RelTol,AbsTol)
%OUTAGEEXACTNUMERICAL Exact numerical evaluation of the one-dimensional outage representation.

    x0 = gammaTh/a;

    F1AtThreshold = aemCDF(x0,mix1);
    L1 = F1AtThreshold^Ns;

    integrand = @(t) transformedIntegrand( ...
        t,x0,a,b,gammaTh,Ns,Nd,mix1,mix2);

    L2 = integral(integrand,0,1, ...
                  'RelTol',RelTol, ...
                  'AbsTol',AbsTol);

    Pout = min(max(L1+L2,0),1);
end

function value = transformedIntegrand( ...
    t,x0,a,b,gammaTh,Ns,Nd,mix1,mix2)
%TRANSFORMEDINTEGRAND Integrand after x=x0+t/(1-t).

    value = zeros(size(t));
    valid = t>0 & t<1;

    if ~any(valid)
        return;
    end

    tValid = t(valid);
    x = x0+tValid./(1-tValid);
    jacobian = 1./(1-tValid).^2;

    denominator = b.*x.*(a.*x-gammaTh);
    yThreshold = gammaTh.*(a.*x+1)./denominator;

    F2 = aemCDF(yThreshold,mix2);
    F1 = aemCDF(x,mix1);
    f1 = aemPDF(x,mix1);

    fX = Ns.*f1.*F1.^(Ns-1);
    FY = F2.^Nd;

    temporary = FY.*fX.*jacobian;
    temporary(~isfinite(temporary)) = 0;

    value(valid) = temporary;
end
