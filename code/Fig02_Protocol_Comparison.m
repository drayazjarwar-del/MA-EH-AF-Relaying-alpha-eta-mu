%% FIGURE 2: TSR, PSR, and IRR outage probability
% Movable-antenna-assisted energy-harvesting AF relaying over
% alpha-eta-mu fading channels.
%
% Analytical curves:
%   Exact numerical evaluation of the one-dimensional outage representation:
%
%   Pout = FX(gamma_th/a)
%        + integral_{gamma_th/a}^{infinity}
%          FY(gamma_th(ax+1)/(bx(ax-gamma_th))) fX(x) dx.
%
% Monte Carlo markers:
%   Parent alpha-eta-mu gains are generated exactly from the Gamma-mixture
%   representation of the parent alpha-eta-mu PDF/CDF.
%
% No Statistics and Machine Learning Toolbox is required; RANDG is used.
%
% MATLAB requirement:
%   Local functions in scripts require MATLAB R2016b or newer.
%
% IMPORTANT THRESHOLD CHOICE
% --------------------------
% cfg.thresholdMode = 'rate':
%   TSR:      gamma_th = 2^(2R/(1-rho)) - 1
%   PSR/IRR:  gamma_th = 2^(2R) - 1
%
% cfg.thresholdMode = 'fixed':
%   All protocols use the same fixed threshold gammaThFixed_dB.
%
% The default is the rate-consistent choice with R = 2 bit/s/Hz.
% To reproduce a common 3-dB threshold comparison,
% change cfg.thresholdMode to 'fixed'.
%
% Generated files:
%   Fig02_Protocol_Comparison.pdf
%   Fig02_Protocol_Comparison.png
%   Fig02_Protocol_Comparison.eps
%   Fig02_Protocol_Comparison_Data.csv

clear; clc; close all;
rng(2026, 'twister');

%% ---------------------------- Configuration -----------------------------
cfg.alpha1 = 2;
cfg.eta1   = 0.5;
cfg.mu1    = 2;
cfg.Omega1 = 1;
cfg.format1 = 1;             % Format I

cfg.alpha2 = 2;
cfg.eta2   = 0.5;
cfg.mu2    = 2;
cfg.Omega2 = 1;
cfg.format2 = 1;             % Format I

cfg.Ns = 3;                  % source candidate positions
cfg.Nd = 3;                  % destination candidate positions

cfg.zeta   = 0.8;            % energy-conversion efficiency
cfg.rho    = 0.3;            % TSR time-switching factor
cfg.lambda = 0.6;            % PSR power-splitting factor

cfg.thresholdMode   = 'rate'; % 'rate' or 'fixed'
cfg.R               = 2;      % target spectral efficiency, bit/s/Hz
cfg.gammaThFixed_dB = 3;      % used only when thresholdMode = 'fixed'

cfg.snr_dB = 0:2:30;          % average source transmit SNR Ps/N0
cfg.Nmc    = 2e8;             % Monte Carlo realizations

% Infinite-series/Gamma-mixture controls
cfg.seriesTol = 1e-13;
cfg.Kmax      = 500;

% Numerical-integration controls
cfg.RelTol = 1e-7;
cfg.AbsTol = 1e-11;

protocols = {'IRR','TSR','PSR'};
lineStyles = {'-','--','-.'};

%% -------------------- Build parent fading distributions -----------------
mix1 = buildAEMMixture(cfg.alpha1, cfg.eta1, cfg.mu1, ...
                       cfg.Omega1, cfg.format1, ...
                       cfg.seriesTol, cfg.Kmax);

mix2 = buildAEMMixture(cfg.alpha2, cfg.eta2, cfg.mu2, ...
                       cfg.Omega2, cfg.format2, ...
                       cfg.seriesTol, cfg.Kmax);

fprintf('First-hop mixture terms retained : %d\n', numel(mix1.w));
fprintf('Second-hop mixture terms retained: %d\n', numel(mix2.w));
fprintf('First-hop retained probability mass : %.15f\n', sum(mix1.w));
fprintf('Second-hop retained probability mass: %.15f\n\n', sum(mix2.w));

%% ---------------------- Basic PDF/CDF validation ------------------------
% These checks should be close to 1. They help detect parameter or
% truncation errors before the outage calculations are performed.
pdfMass1 = integral(@(z) aemPDF(z, mix1), 0, Inf, ...
                    'RelTol', 1e-8, 'AbsTol', 1e-11);
pdfMass2 = integral(@(z) aemPDF(z, mix2), 0, Inf, ...
                    'RelTol', 1e-8, 'AbsTol', 1e-11);

fprintf('Integral of first-hop parent PDF : %.12f\n', pdfMass1);
fprintf('Integral of second-hop parent PDF: %.12f\n\n', pdfMass2);

if abs(pdfMass1 - 1) > 1e-6 || abs(pdfMass2 - 1) > 1e-6
    warning(['A parent PDF does not integrate sufficiently close to one. ', ...
             'Increase cfg.Kmax or inspect the fading parameters.']);
end

%% ---------------- Generate selected channel gains once -----------------
% The same independent fading realizations are reused for every SNR and
% protocol. This reduces Monte Carlo noise in protocol comparisons.

fprintf('Generating %g selected first-hop gains...\n', cfg.Nmc);
X = zeros(cfg.Nmc,1);
for n = 1:cfg.Ns
    X = max(X, sampleAEM(mix1, cfg.Nmc));
end

fprintf('Generating %g selected second-hop gains...\n\n', cfg.Nmc);
Y = zeros(cfg.Nmc,1);
for m = 1:cfg.Nd
    Y = max(Y, sampleAEM(mix2, cfg.Nmc));
end

%% ---------------- Analytical and Monte Carlo outage --------------------
snrLin = 10.^(cfg.snr_dB/10);
nSNR = numel(snrLin);
nProt = numel(protocols);

PoutAna = zeros(nSNR, nProt);
PoutMC  = zeros(nSNR, nProt);
outageCounts = zeros(nSNR, nProt);

for p = 1:nProt
    protocol = protocols{p};
    fprintf('Evaluating %s...\n', protocol);

    for q = 1:nSNR
        [a, b, gammaTh] = protocolParameters( ...
            snrLin(q), protocol, cfg);

        % Exact numerical evaluation of the one-dimensional outage representation
        PoutAna(q,p) = outageExactNumerical( ...
            a, b, gammaTh, cfg.Ns, cfg.Nd, ...
            mix1, mix2, cfg.RelTol, cfg.AbsTol);

        % Independent Monte Carlo evaluation of Eq. (13)
        gamma1  = a .* X;
        gamma2  = b .* X .* Y;
        gammaAF = (gamma1 .* gamma2) ./ (gamma1 + gamma2 + 1);

        outageCounts(q,p) = sum(gammaAF < gammaTh);
        PoutMC(q,p) = outageCounts(q,p) / cfg.Nmc;

        fprintf(['  SNR = %5.1f dB | gamma_th = %10.4g | ', ...
                 'Analytical = %.4e | MC = %.4e | outages = %d\n'], ...
                 cfg.snr_dB(q), gammaTh, PoutAna(q,p), ...
                 PoutMC(q,p), outageCounts(q,p));
    end
    fprintf('\n');
end

%% ------------------------------ Plot -----------------------------------
fig = figure('Color','w', ...
    'Units','inches', ...
    'Position',[1 1 3.5 3.10], ...
    'Renderer','opengl');

ax = axes(fig);
set(ax,'YScale','log');
hold(ax,'on');

hLine = gobjects(nProt,1);

lineColors = [
    0       0.4470 0.7410;   % IRR - blue
    0.8500 0.3250 0.0980;    % TSR - orange/red
    0.4660 0.6740 0.1880     % PSR - green
];

for p = 1:nProt

    % Analytical curves
    hLine(p) = plot(ax,cfg.snr_dB,PoutAna(:,p), ...
        lineStyles{p}, ...
        'Color',lineColors(p,:), ...
        'LineWidth',1, ...
        'DisplayName',protocols{p});

    % Simulation markers
    mcToPlot = PoutMC(:,p);
    mcToPlot(mcToPlot==0) = NaN;

    plot(ax,cfg.snr_dB,mcToPlot,'o', ...
        'LineStyle','none', ...
        'Color','k', ...
        'MarkerFaceColor','none', ...
        'MarkerSize',5, ...
        'LineWidth',0.8, ...
        'HandleVisibility','off');
end

% Simulation legend marker
hSim = plot(ax,NaN,NaN,'ko', ...
    'LineStyle','none', ...
    'MarkerFaceColor','none', ...
    'MarkerSize',5, ...
    'LineWidth',0.8, ...
    'DisplayName','Simulation');

grid(ax,'on');
box(ax,'on');

ax.YMinorGrid = 'on';
ax.XMinorGrid = 'off';
ax.GridAlpha = 0.22;
ax.MinorGridAlpha = 0.12;

ax.FontName = 'Times New Roman';
ax.FontSize = 8;
ax.LineWidth = 0.8;
ax.TickDir = 'in';
ax.Layer = 'top';

xlabel(ax,'Average source transmit SNR', ...
    'FontSize',9);

ylabel(ax,'Outage probability', ...
    'FontSize',9);

lgd = legend(ax,[hLine;hSim], ...
    {'IRR','TSR','PSR','Simulation'}, ...
    'Location','southwest', ...
    'Interpreter','latex');

lgd.FontName = 'Times New Roman';
lgd.FontSize = 11;
lgd.Box = 'on';
lgd.LineWidth = 1;
lgd.ItemTokenSize = [11 11];

xlim(ax,[min(cfg.snr_dB) max(cfg.snr_dB)]);
ylim(ax,[1e-6 1]);

%% ---------------------------- Save outputs ------------------------------
exportgraphics(fig, 'Fig02_Protocol_Comparison.pdf', ...
               'ContentType','vector', ...
               'BackgroundColor','white');
exportgraphics(fig, 'Fig02_Protocol_Comparison.png', ...
               'Resolution',600);
print(fig, 'Fig02_Protocol_Comparison.eps', '-depsc2', '-r600');

% Save numerical results for later verification and manuscript tables.
T = table(cfg.snr_dB(:), ...
          PoutAna(:,1), PoutMC(:,1), outageCounts(:,1), ...
          PoutAna(:,2), PoutMC(:,2), outageCounts(:,2), ...
          PoutAna(:,3), PoutMC(:,3), outageCounts(:,3), ...
    'VariableNames', {'SNR_dB', ...
                      'IRR_Analytical','IRR_MC','IRR_OutageCount', ...
                      'TSR_Analytical','TSR_MC','TSR_OutageCount', ...
                      'PSR_Analytical','PSR_MC','PSR_OutageCount'});
writetable(T, 'Fig02_Protocol_Comparison_Data.csv');

fprintf('Figure and numerical results have been exported.\n');

%% =======================================================================
%                           LOCAL FUNCTIONS
% ========================================================================

function mix = buildAEMMixture(alpha, eta, mu, Omega, formatNo, tol, Kmax)
%BUILDAEMMIXTURE Build the Gamma-mixture representation of alpha-eta-mu.
%
% For the parent gain Z, define
%   U = c (Z/Omega)^(alpha/2),  c = 2 mu h.
%
% The series PDF/CDF can be written as a positive mixture:
%   K ~ categorical(w_k)
%   U | K=k ~ Gamma(2(mu+k), 1)
%   Z = Omega (U/c)^(2/alpha).
%
% The CDF is therefore
%   F_Z(z) = sum_k w_k P(2(mu+k), c(z/Omega)^(alpha/2)),
% where P(.,.) is the regularized lower incomplete gamma function.

    validateattributes(alpha, {'numeric'}, {'scalar','real','positive'});
    validateattributes(mu,    {'numeric'}, {'scalar','real','positive'});
    validateattributes(Omega, {'numeric'}, {'scalar','real','positive'});

    [h,H] = etaAuxiliary(eta, formatNo);
    c = 2*mu*h;

    if abs(H) < 1e-15
        k = 0;
        w = 1;
    else
        k = (0:Kmax).';

        logw = log(2) + 0.5*log(pi) ...
             + 2*mu*log(mu) + mu*log(h) - gammaln(mu) ...
             + 2*k.*log(abs(mu*H)) ...
             - gammaln(k+1) ...
             - gammaln(mu+k+0.5) ...
             - 2*(mu+k).*log(c) ...
             + gammaln(2*(mu+k));

        wAll = exp(logw);
        totalMass = sum(wAll);

        if abs(totalMass - 1) > 1e-8
            warning(['The series mass at Kmax=%d is %.12g rather than 1. ', ...
                     'Increase Kmax for these fading parameters.'], ...
                     Kmax, totalMass);
        end

        cumulativeMass = cumsum(wAll);
        last = find(cumulativeMass >= 1-tol, 1, 'first');

        if isempty(last)
            last = numel(wAll);
        end

        k = k(1:last);
        w = wAll(1:last);

        % Renormalize only the numerically negligible truncated tail.
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

function [h,H] = etaAuxiliary(eta, formatNo)
%ETAAUXILIARY Auxiliary h and H parameters for the eta-mu formats.

    switch formatNo
        case 1
            if ~(eta > 0)
                error('For Format I, eta must satisfy eta > 0.');
            end
            h = (1+eta)^2/(4*eta);
            H = (1-eta^2)/(4*eta);

        case 2
            if ~(eta > -1 && eta < 1)
                error('For Format II, eta must satisfy -1 < eta < 1.');
            end
            h = 1/(1-eta^2);
            H = eta/(1-eta^2);

        otherwise
            error('formatNo must be 1 or 2.');
    end
end

function F = aemCDF(z, mix)
%AEMCDF Parent alpha-eta-mu power-gain CDF.

    originalSize = size(z);
    z = z(:).';

    F = zeros(size(z));
    positive = z > 0 & isfinite(z);

    if any(positive)
        u = mix.c .* (z(positive)./mix.Omega).^(mix.alpha/2);
        Fpos = zeros(size(u));

        for k = 1:numel(mix.w)
            Fpos = Fpos + mix.w(k) .* ...
                gammainc(u, mix.shape(k), 'lower');
        end

        F(positive) = Fpos;
    end

    F(isinf(z) & z > 0) = 1;
    F(z <= 0) = 0;
    F = min(max(F,0),1);
    F = reshape(F, originalSize);
end

function f = aemPDF(z, mix)
%AEMPDF Parent alpha-eta-mu power-gain PDF.

    originalSize = size(z);
    z = z(:).';

    f = zeros(size(z));
    positive = z > 0 & isfinite(z);

    if any(positive)
        zp = z(positive);
        u = mix.c .* (zp./mix.Omega).^(mix.alpha/2);

        duDz = mix.c .* (mix.alpha/2) ./ mix.Omega ...
             .* (zp./mix.Omega).^(mix.alpha/2 - 1);

        fpos = zeros(size(zp));

        for k = 1:numel(mix.w)
            s = mix.shape(k);
            gammaPDF = exp((s-1).*log(u) - u - gammaln(s));
            fpos = fpos + mix.w(k).*gammaPDF.*duDz;
        end

        f(positive) = fpos;
    end

    % For the present manuscript parameters the density tends to zero at
    % z=0. The exact point value does not affect numerical integration.
    f(z <= 0) = 0;
    f = reshape(f, originalSize);
end

function z = sampleAEM(mix, N)
%SAMPLEAEM Generate N exact alpha-eta-mu parent-gain samples.

    uCat = rand(N,1);
    cumulativeWeights = cumsum(mix.w);
    cumulativeWeights(end) = 1;

    % Categorical mixture index without requiring RANDSAMPLE.
    idx = ones(N,1,'uint16');
    for k = 1:numel(cumulativeWeights)-1
        idx(uCat > cumulativeWeights(k)) = k+1;
    end

    shapePerSample = mix.shape(double(idx));

    % RANDG generates Gamma(shape, scale=1) random variables and is part of
    % base MATLAB.
    gammaSample = randg(shapePerSample);

    z = mix.Omega .* (gammaSample./mix.c).^(2/mix.alpha);
end

function [a,b,gammaTh] = protocolParameters(snrLin, protocol, cfg)
%PROTOCOLPARAMETERS Unified a, b, and outage-threshold parameters.

    switch upper(protocol)
        case 'TSR'
            a = snrLin;
            b = (2*cfg.zeta*cfg.rho/(1-cfg.rho))*snrLin;

        case 'PSR'
            a = (1-cfg.lambda)*snrLin;
            b = cfg.zeta*cfg.lambda*snrLin;

        case 'IRR'
            a = snrLin;
            b = cfg.zeta*snrLin;

        otherwise
            error('Unknown protocol: %s', protocol);
    end

    switch lower(cfg.thresholdMode)
        case 'rate'
            if strcmpi(protocol,'TSR')
                gammaTh = 2^(2*cfg.R/(1-cfg.rho)) - 1;
            else
                gammaTh = 2^(2*cfg.R) - 1;
            end

        case 'fixed'
            gammaTh = 10^(cfg.gammaThFixed_dB/10);

        otherwise
            error('cfg.thresholdMode must be ''rate'' or ''fixed''.');
    end
end

function Pout = outageExactNumerical(a, b, gammaTh, Ns, Nd, ...
                                     mix1, mix2, RelTol, AbsTol)
%OUTAGEEXACTNUMERICAL Exact numerical evaluation of the one-dimensional outage representation.
%
% The transformation
%   x = x0 + t/(1-t),  0 < t < 1,
% maps [x0,infinity) to [0,1) and avoids evaluating directly at infinity.

    x0 = gammaTh/a;

    F1x0 = aemCDF(x0, mix1);
    L1 = F1x0^Ns;

    integrand = @(t) transformedIntegrand( ...
        t, x0, a, b, gammaTh, Ns, Nd, mix1, mix2);

    L2 = integral(integrand, 0, 1, ...
                  'RelTol', RelTol, ...
                  'AbsTol', AbsTol);

    Pout = min(max(L1 + L2, 0), 1);
end

function val = transformedIntegrand(t, x0, a, b, gammaTh, ...
                                    Ns, Nd, mix1, mix2)
%TRANSFORMEDINTEGRAND Vectorized integrand after x=x0+t/(1-t).

    val = zeros(size(t));
    valid = t > 0 & t < 1;

    if ~any(valid)
        return;
    end

    tv = t(valid);
    x = x0 + tv./(1-tv);
    jacobian = 1./(1-tv).^2;

    denominator = b.*x.*(a.*x-gammaTh);
    yThreshold = gammaTh.*(a.*x+1)./denominator;

    F2 = aemCDF(yThreshold, mix2);
    F1 = aemCDF(x, mix1);
    f1 = aemPDF(x, mix1);

    fX = Ns.*f1.*F1.^(Ns-1);
    FY = F2.^Nd;

    temp = FY.*fX.*jacobian;
    temp(~isfinite(temp)) = 0;

    val(valid) = temp;
end
