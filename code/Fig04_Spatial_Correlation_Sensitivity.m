%% FIGURE 4: Spatial-correlation sensitivity of MA selection
% MA-assisted EH variable-gain AF relaying over alpha-eta-mu fading.
%
% Purpose
% -------
% Quantify the degradation of the i.i.d. candidate-position benchmark when
% candidate MA positions are spatially correlated. The correlation-aware
% Monte Carlo model preserves the marginal Format-I alpha-eta-mu law at
% every candidate position (for 2*mu integer) while imposing Clarke's
% isotropic 2-D spatial correlation on the underlying scattered Gaussian
% components:
%
%   R[p,q] = J_0(2*pi*|p-q|*Delta/lambda_c).
%
% With N candidate positions, the movement aperture is
%
%   W = (N-1)*Delta.
%
% The i.i.d. curve is evaluated analytically from the exact 1-D outage
% expression used in the submitted manuscript. Correlated curves are Monte
% Carlo results. A fixed-antenna analytical benchmark (Ns=Nd=1) is included.
%
% Default sensitivity setting
% ---------------------------
%   alpha1=alpha2=2, eta1=eta2=0.5, mu1=mu2=2, Omega1=Omega2=1
%   PSR: lambda=0.6, zeta=0.8, R=2 bit/s/Hz
%   Ns=Nd=3
%   Delta/lambda_c = 0.05, 0.10, 0.20
%   target outage for SNR-penalty table = 1e-3
%
% The correlated sample generator uses no Statistics and Machine Learning
% Toolbox functions.
%
% Outputs
% -------
%   Fig04_Spatial_Correlation_Sensitivity.pdf
%   Fig04_Spatial_Correlation_Sensitivity.png
%   Fig04_Spatial_Correlation_Sensitivity.csv
%
% MATLAB requirement: local functions in scripts require R2016b or newer.

clear; clc; close all;
rng(20260914,'twister');

%% ---------------------------- Configuration -----------------------------
cfg.alpha1 = 2;
cfg.eta1   = 0.5;
cfg.mu1    = 2;
cfg.Omega1 = 1;
cfg.format1 = 1;                 % correlation generator below: Format I

cfg.alpha2 = 2;
cfg.eta2   = 0.5;
cfg.mu2    = 2;
cfg.Omega2 = 1;
cfg.format2 = 1;

cfg.Ns = 3;
cfg.Nd = 3;

cfg.zeta   = 0.8;
cfg.lambda = 0.6;
cfg.R      = 2;                  % bit/s/Hz
cfg.gammaTh = 2^(2*cfg.R)-1;

cfg.snr_dB = 15:1:30;
cfg.targetOutage = 1e-3;

% Candidate spacing normalized by carrier wavelength lambda_c.
cfg.deltaNorm = [0.05 0.10 0.20];

% Monte Carlo settings. 2e7 gives about 2e4 tail samples at Pout=1e-3.
% For a quick test, temporarily use cfg.Nmc = 2e5.
cfg.Nmc = 2e7;
cfg.chunkSize = 1e5;

% Streaming histogram for the 99.9th percentile of required SNR.
cfg.histMin_dB  = 5;
cfg.histMax_dB  = 50;
cfg.histStep_dB = 0.005;
cfg.histEdges_dB = cfg.histMin_dB:cfg.histStep_dB:cfg.histMax_dB;

% Exact i.i.d. analytical settings.
cfg.seriesTol = 1e-13;
cfg.Kmax      = 500;
cfg.RelTol    = 1e-8;
cfg.AbsTol    = 1e-12;
cfg.snrTol_dB = 1e-4;

%% -------------------- Build exact marginal distributions ----------------
mix1 = buildAEMMixture(cfg.alpha1,cfg.eta1,cfg.mu1,cfg.Omega1, ...
                       cfg.format1,cfg.seriesTol,cfg.Kmax);
mix2 = buildAEMMixture(cfg.alpha2,cfg.eta2,cfg.mu2,cfg.Omega2, ...
                       cfg.format2,cfg.seriesTol,cfg.Kmax);

fprintf('\nSpatial-correlation sensitivity study\n');
fprintf('  Ns = Nd = %d, lambda = %.2f, target outage = %.1e\n', ...
        cfg.Ns,cfg.lambda,cfg.targetOutage);
fprintf('  Monte Carlo realizations per spacing = %.0e\n\n',cfg.Nmc);

%% ------------------- Exact i.i.d. and fixed benchmarks ------------------
nSNR = numel(cfg.snr_dB);
PoutIID = zeros(1,nSNR);
PoutFixed = zeros(1,nSNR);

for k = 1:nSNR
    snrLin = 10^(cfg.snr_dB(k)/10);
    a = (1-cfg.lambda)*snrLin;
    b = cfg.zeta*cfg.lambda*snrLin;

    PoutIID(k) = outageExactNumerical( ...
        a,b,cfg.gammaTh,cfg.Ns,cfg.Nd,mix1,mix2,cfg.RelTol,cfg.AbsTol);

    PoutFixed(k) = outageExactNumerical( ...
        a,b,cfg.gammaTh,1,1,mix1,mix2,cfg.RelTol,cfg.AbsTol);
end

snrIID_target_dB = requiredSNRExact(cfg.targetOutage,cfg.Ns,cfg.Nd,cfg,mix1,mix2);
snrFixed_target_dB = requiredSNRExact(cfg.targetOutage,1,1,cfg,mix1,mix2);

%% ---------------- Correlation-aware Monte Carlo simulation --------------
nDelta = numel(cfg.deltaNorm);
PoutCorr = zeros(nDelta,nSNR);
requiredSNRCorr_dB = zeros(nDelta,1);
nearestCorr = zeros(nDelta,1);
nearestGainCorr = zeros(nDelta,1);
apertureNorm = zeros(nDelta,1);

nBins = numel(cfg.histEdges_dB)-1;

for dIdx = 1:nDelta
    deltaNorm = cfg.deltaNorm(dIdx);

    RspS = clarkeCorrelationMatrix(cfg.Ns,deltaNorm);
    RspD = clarkeCorrelationMatrix(cfg.Nd,deltaNorm);
    nearestCorr(dIdx) = RspS(1,min(2,cfg.Ns));
    % For the default alpha=2 Format-I Gaussian-cluster construction,
    % the Pearson correlation of the candidate power gains is r_field^2.
    nearestGainCorr(dIdx) = nearestCorr(dIdx)^2;
    apertureNorm(dIdx) = (cfg.Ns-1)*deltaNorm;

    fprintf(['Delta/lambda_c = %.3f | W/lambda_c = %.3f | ', ...
             'nearest field corr = %.6f | nearest gain corr = %.6f\n'], ...
            deltaNorm,apertureNorm(dIdx),nearestCorr(dIdx),nearestGainCorr(dIdx));

    outageCounts = zeros(1,nSNR);
    histCounts = zeros(nBins,1);
    underflowCount = 0;
    overflowCount = 0;
    processed = 0;

    while processed < cfg.Nmc
        thisChunk = min(cfg.chunkSize,cfg.Nmc-processed);

        X = sampleCorrelatedSelectedAEMFormatI( ...
            thisChunk,cfg.Ns,RspS,cfg.alpha1,cfg.eta1,cfg.mu1,cfg.Omega1);
        Y = sampleCorrelatedSelectedAEMFormatI( ...
            thisChunk,cfg.Nd,RspD,cfg.alpha2,cfg.eta2,cfg.mu2,cfg.Omega2);

        sReq_dB = requiredSNRPerRealizationPSR( ...
            X,Y,cfg.lambda,cfg.zeta,cfg.gammaTh);

        for k = 1:nSNR
            outageCounts(k) = outageCounts(k) + sum(sReq_dB > cfg.snr_dB(k));
        end

        underflowCount = underflowCount + sum(sReq_dB < cfg.histMin_dB);
        overflowCount  = overflowCount  + sum(sReq_dB >= cfg.histMax_dB);
        histCounts = histCounts + histcounts(sReq_dB,cfg.histEdges_dB).';

        processed = processed + thisChunk;

        if mod(processed,1e6)==0 || processed==cfg.Nmc
            fprintf('  processed %.0f / %.0f\n',processed,cfg.Nmc);
        end
    end

    PoutCorr(dIdx,:) = outageCounts/cfg.Nmc;
    requiredSNRCorr_dB(dIdx) = quantileFromHistogram( ...
        histCounts,underflowCount,overflowCount,cfg.histEdges_dB, ...
        cfg.Nmc,1-cfg.targetOutage);

    fprintf('  required SNR at Pout=1e-3: %.4f dB\n\n', ...
            requiredSNRCorr_dB(dIdx));
end

%% ----------------------------- Diagnostics ------------------------------
penaltyVsIID_dB = requiredSNRCorr_dB-snrIID_target_dB;

fprintf('============================================================\n');
fprintf('Required SNR at Pout = %.1e\n',cfg.targetOutage);
fprintf('  Fixed antenna (Ns=Nd=1): %.4f dB\n',snrFixed_target_dB);
fprintf('  i.i.d. MA benchmark:      %.4f dB\n',snrIID_target_dB);
for dIdx = 1:nDelta
    fprintf(['  Delta/lambda_c=%.2f, W/lambda_c=%.2f, r1=%.4f: ', ...
             'gain corr=%.4f: %.4f dB (penalty %.4f dB vs i.i.d.)\n'], ...
             cfg.deltaNorm(dIdx),apertureNorm(dIdx),nearestCorr(dIdx),nearestGainCorr(dIdx), ...
             requiredSNRCorr_dB(dIdx),penaltyVsIID_dB(dIdx));
end
fprintf('============================================================\n\n');

%% ----------------------------- Save data --------------------------------
% One row per SNR for the outage curves.
Tcurve = table(cfg.snr_dB(:),PoutIID(:),PoutFixed(:), ...
    PoutCorr(1,:).',PoutCorr(2,:).',PoutCorr(3,:).', ...
    'VariableNames',{'SNR_dB','Pout_iid_exact','Pout_fixed_exact', ...
    'Pout_corr_Delta005','Pout_corr_Delta010','Pout_corr_Delta020'});

writetable(Tcurve,'Fig04_Spatial_Correlation_Sensitivity.csv');

Tsummary = table(cfg.deltaNorm(:),apertureNorm,nearestCorr,nearestGainCorr, ...
    requiredSNRCorr_dB,penaltyVsIID_dB, ...
    'VariableNames',{'Delta_over_lambda','W_over_lambda','NearestFieldCorrelation', ...
    'NearestGainCorrelation','RequiredSNR_dB','PenaltyVsIID_dB'});

disp(Tsummary);

%% ------------------------ Publication-ready figure ----------------------
figure('Color','w','Position',[120 80 760 580]);
hold on;

% Exact benchmarks.
semilogy(cfg.snr_dB,PoutFixed,':','LineWidth',1.8, ...
    'DisplayName','$N_s=N_d=1$ (exact)');
semilogy(cfg.snr_dB,PoutIID,'-','LineWidth',2.2, ...
    'DisplayName','$N_s=N_d=3$, i.i.d. (exact)');

styles = {'--','-.','-'};
for dIdx = 1:nDelta
    valid = PoutCorr(dIdx,:)>0;
    semilogy(cfg.snr_dB(valid),PoutCorr(dIdx,valid),styles{dIdx}, ...
        'LineWidth',1.6, ...
        'Marker','o','MarkerSize',5.5, ...
        'MarkerFaceColor','none', ...
        'DisplayName',sprintf('$\\Delta/\\lambda_c=%.2f$ (Simulation)', ...
                              cfg.deltaNorm(dIdx)));
end

xlabel('Average source transmit SNR (dB)','Interpreter','latex');
ylabel('Outage probability','Interpreter','latex');
xlim([min(cfg.snr_dB) max(cfg.snr_dB)]);
ylim([1e-5 1]);
grid on; box on;
set(gca,'YScale','log','FontName','Times New Roman','FontSize',12, ...
        'LineWidth',1.0);
legend('Location','southwest','Interpreter','latex','FontSize',10);

exportgraphics(gcf,'Fig04_Spatial_Correlation_Sensitivity.pdf', ...
               'ContentType','vector');
exportgraphics(gcf,'Fig04_Spatial_Correlation_Sensitivity.png', ...
               'Resolution',600);

%% ========================================================================
%                              LOCAL FUNCTIONS
% ========================================================================

function R = clarkeCorrelationMatrix(N,deltaNorm)
%CLARKECORRELATIONMATRIX Clarke/Jakes spatial field correlation.
% deltaNorm = Delta/lambda_c.

    idx = 0:N-1;
    D = abs(idx(:)-idx(:).');
    R = besselj(0,2*pi*deltaNorm*D);
    R = 0.5*(R+R.');

    % Numerical safeguard only; the selected spacings give positive-
    % definite matrices for the dimensions used here.
    minEig = min(eig(R));
    if minEig < -1e-10
        error('Spatial correlation matrix is not positive semidefinite.');
    elseif minEig < 1e-12
        R = R + (1e-12-minEig)*eye(N);
    end
end

function Zsel = sampleCorrelatedSelectedAEMFormatI( ...
    N,Npos,Rsp,alpha,eta,mu,Omega)
%SAMPLECORRELATEDSELECTEDAEMFORMATI Correlated Format-I alpha-eta-mu gains.
%
% The construction is exact at each candidate position when 2*mu is an
% integer. The underlying in-phase and quadrature Gaussian components are
% spatially correlated according to Rsp. Their variances satisfy the
% Format-I power ratio eta. The dimensionless eta-mu power W has E[W]=1;
% the alpha-nonlinearity transformation is Z = Omega*W^(2/alpha).

    nu = round(2*mu);
    if abs(nu-2*mu)>1e-12 || nu<1
        error(['Correlation-aware physical generator requires 2*mu to ', ...
               'be a positive integer. The paper default mu=2 satisfies this.']);
    end
    if eta<=0
        error('Format-I eta must be positive.');
    end

    L = chol(Rsp,'lower');

    varI = eta/(nu*(1+eta));
    varQ = 1/(nu*(1+eta));

    W = zeros(N,Npos);
    for r = 1:nu
        GI = sqrt(varI)*(randn(N,Npos)*L.');
        GQ = sqrt(varQ)*(randn(N,Npos)*L.');
        W = W + GI.^2 + GQ.^2;
    end

    Zcand = Omega.*W.^(2/alpha);
    Zsel = max(Zcand,[],2);
end

function sReq_dB = requiredSNRPerRealizationPSR(X,Y,lambda,zeta,gammaTh)
%REQUIREDSNRPERREALIZATIONPSR Positive root of gamma_AF=gamma_th.

    A = 1-lambda;
    B = zeta*lambda;

    quadraticCoefficient = A*B.*X.^2.*Y;
    linearMagnitude = gammaTh.*X.*(A+B.*Y);

    sReq = (linearMagnitude + ...
        sqrt(linearMagnitude.^2 + 4.*quadraticCoefficient.*gammaTh)) ./ ...
        (2.*quadraticCoefficient);

    sReq_dB = 10.*log10(sReq);
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
            error('Requested quantile lies above histogram maximum.');
        else
            error('Unable to determine Monte Carlo quantile.');
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

    fractionWithinBin = (targetRank-countBefore)/countInBin;
    fractionWithinBin = min(max(fractionWithinBin,0),1);
    q_dB = edges_dB(binIndex) + fractionWithinBin* ...
           (edges_dB(binIndex+1)-edges_dB(binIndex));
end

function snrReq_dB = requiredSNRExact(targetOutage,Ns,Nd,cfg,mix1,mix2)
%REQUIREDSNREXACT Bisection on the exact i.i.d. outage expression.

    lo = 0;
    hi = 45;

    fLo = exactPoutAtSNR(lo,Ns,Nd,cfg,mix1,mix2)-targetOutage;
    fHi = exactPoutAtSNR(hi,Ns,Nd,cfg,mix1,mix2)-targetOutage;

    if fLo<0 || fHi>0
        error('Bisection bracket does not contain the target outage.');
    end

    while (hi-lo)>cfg.snrTol_dB
        mid = 0.5*(lo+hi);
        fMid = exactPoutAtSNR(mid,Ns,Nd,cfg,mix1,mix2)-targetOutage;
        if fMid>0
            lo = mid;
        else
            hi = mid;
        end
    end
    snrReq_dB = 0.5*(lo+hi);
end

function Pout = exactPoutAtSNR(snr_dB,Ns,Nd,cfg,mix1,mix2)
    s = 10^(snr_dB/10);
    a = (1-cfg.lambda)*s;
    b = cfg.zeta*cfg.lambda*s;
    Pout = outageExactNumerical( ...
        a,b,cfg.gammaTh,Ns,Nd,mix1,mix2,cfg.RelTol,cfg.AbsTol);
end

function Pout = outageExactNumerical(a,b,gammaTh,Ns,Nd,mix1,mix2,RelTol,AbsTol)
%OUTAGEEXACTNUMERICAL Exact one-dimensional outage representation.

    x0 = gammaTh/a;
    F1x0 = aemCDF(x0,mix1);
    L1 = F1x0^Ns;

    integrand = @(t) transformedIntegrand( ...
        t,x0,a,b,gammaTh,Ns,Nd,mix1,mix2);

    L2 = integral(integrand,0,1,'RelTol',RelTol,'AbsTol',AbsTol);
    Pout = min(max(L1+L2,0),1);
end

function val = transformedIntegrand(t,x0,a,b,gammaTh,Ns,Nd,mix1,mix2)
    val = zeros(size(t));
    inside = t>=0 & t<1;
    if ~any(inside)
        return;
    end

    ti = t(inside);
    x = x0 + ti./(1-ti);

    yThreshold = gammaTh.*(a.*x+1) ./ ...
        (b.*x.*(a.*x-gammaTh));

    F2 = aemCDF(yThreshold,mix2);
    F1 = aemCDF(x,mix1);
    f1 = aemPDF(x,mix1);

    fX = Ns.*F1.^(Ns-1).*f1;
    G = F2.^Nd;

    val(inside) = G.*fX./(1-ti).^2;
end

function mix = buildAEMMixture(alpha,eta,mu,Omega,formatNo,tol,Kmax)
%BUILDAEMMIXTURE Gamma-mixture representation of alpha-eta-mu fading.

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
    originalSize = size(z);
    z = z(:).';
    F = zeros(size(z));
    positive = z>0 & isfinite(z);

    if any(positive)
        u = mix.c.*(z(positive)./mix.Omega).^(mix.alpha/2);
        Fpos = zeros(size(u));
        for k = 1:numel(mix.w)
            Fpos = Fpos + mix.w(k).*gammainc(u,mix.shape(k),'lower');
        end
        F(positive) = Fpos;
    end

    F(isinf(z) & z>0) = 1;
    F(z<=0) = 0;
    F = min(max(F,0),1);
    F = reshape(F,originalSize);
end

function f = aemPDF(z,mix)
    originalSize = size(z);
    z = z(:).';
    f = zeros(size(z));
    positive = z>0 & isfinite(z);

    if any(positive)
        zp = z(positive);
        u = mix.c.*(zp./mix.Omega).^(mix.alpha/2);
        duDz = mix.c.*(mix.alpha/2)./mix.Omega .* ...
               (zp./mix.Omega).^(mix.alpha/2-1);
        fpos = zeros(size(zp));

        for k = 1:numel(mix.w)
            s = mix.shape(k);
            gammaPDF = exp((s-1).*log(u)-u-gammaln(s));
            fpos = fpos + mix.w(k).*gammaPDF.*duDz;
        end
        f(positive) = fpos;
    end

    f(z<=0) = 0;
    f = reshape(f,originalSize);
end
