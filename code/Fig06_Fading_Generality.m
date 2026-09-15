%% FIGURE 6: Generality of the alpha-eta-mu framework
% Reproduces Fig. 6 of the revised manuscript.
% This single figure demonstrates:
%   (i) the default homogeneous alpha-eta-mu case,
%   (ii) sensitivity to alpha,
%   (iii) sensitivity to eta, and
%   (iv) a genuinely heterogeneous Rayleigh/Nakagami-m two-hop case.
%
% All curves use PSR with Ns = Nd = 3, lambda = 0.6, zeta = 0.8, R = 2.
% Analytical curves use the exact one-dimensional outage representation.
% Monte Carlo markers are shown only where >=20 outage events are observed.
% For the Rayleigh/Nakagami-m case, samples are generated directly from
% exponential/Gamma power-gain laws, independently of the general sampler.
%
% Generated files:
%   Fig06_Fading_Generality.pdf   (direct replacement in the LaTeX project)
%   Fig06_Fading_Generality.png
%   Fig06_Fading_Generality_Data.csv

clear; clc; close all;
rng(2026,'twister');

%% ---------------------------- Configuration -----------------------------
cfg.zeta   = 0.8;
cfg.lambda = 0.6;
cfg.R      = 2;
cfg.Ns     = 3;
cfg.Nd     = 3;
cfg.Omega1 = 1;
cfg.Omega2 = 1;
cfg.format1 = 1;   % Format I
cfg.format2 = 1;

cfg.snr_dB = 0:2:30;
cfg.Nmc = 1e8;          % very-low-outage validation
cfg.chunkSize = 1e6;    % memory-efficient chunking
cfg.minEventsToPlot = 20;

cfg.seriesTol = 1e-13;
cfg.Kmax = 500;
cfg.RelTol = 1e-7;
cfg.AbsTol = 1e-11;

% Columns: alpha1 eta1 mu1 alpha2 eta2 mu2
% Case 1: default homogeneous alpha-eta-mu.
% Case 2: alpha sensitivity (both hops, alpha reduced from 2 to 1.5).
% Case 3: eta sensitivity (both hops, eta reduced from 0.5 to 0.2).
% Case 4: genuinely heterogeneous Rayleigh/Nakagami-m.
%         Rayleigh: alpha=2, eta=1, mu=0.5.
%         Nakagami-m with m=2: alpha=2, eta=1, mu=m/2=1.
cases = [ ...
    2.0 0.5 2.0   2.0 0.5 2.0; ...
    1.5 0.5 2.0   1.5 0.5 2.0; ...
    2.0 0.2 2.0   2.0 0.2 2.0; ...
    2.0 1.0 0.5   2.0 1.0 1.0  ...
    ];

caseLabels = { ...
    'Default: $\alpha=2,\eta=0.5,\mu=2$', ...
    '$\alpha=1.5$ (homogeneous)', ...
    '$\eta=0.2$ (homogeneous)', ...
    'Rayleigh/Nakagami-$m$ ($m=2$)' ...
    };

lineStyles = {'-','--','-.',':'};
lineColors = [ ...
    0.0000 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.4660 0.6740 0.1880; ...
    0.4940 0.1840 0.5560  ...
    ];

%% ---------------------- PSR protocol parameters -------------------------
snrLin  = 10.^(cfg.snr_dB/10);
gammaTh = 2^(2*cfg.R)-1;
aVec = (1-cfg.lambda).*snrLin;
bVec = cfg.zeta*cfg.lambda.*snrLin;

nSNR = numel(snrLin);
nCases = size(cases,1);
PoutAna = zeros(nSNR,nCases);
PoutMC = zeros(nSNR,nCases);
outageCounts = zeros(nSNR,nCases);

%% ---------------- Build distributions + analytical curves --------------
mix1 = cell(nCases,1);
mix2 = cell(nCases,1);

fprintf('Figure 6: heterogeneous/generalised fading study\n');
fprintf('Building parent distributions and exact outage curves...\n');

for c = 1:nCases
    alpha1 = cases(c,1); eta1 = cases(c,2); mu1 = cases(c,3);
    alpha2 = cases(c,4); eta2 = cases(c,5); mu2 = cases(c,6);

    mix1{c} = buildAEMMixture(alpha1,eta1,mu1,cfg.Omega1,cfg.format1, ...
        cfg.seriesTol,cfg.Kmax);
    mix2{c} = buildAEMMixture(alpha2,eta2,mu2,cfg.Omega2,cfg.format2, ...
        cfg.seriesTol,cfg.Kmax);

    fprintf('\nCase %d: (a1,e1,m1)=(%.2f,%.2f,%.2f), ', ...
        c,alpha1,eta1,mu1);
    fprintf('(a2,e2,m2)=(%.2f,%.2f,%.2f)\n',alpha2,eta2,mu2);

    for q = 1:nSNR
        PoutAna(q,c) = outageExactNumerical( ...
            aVec(q),bVec(q),gammaTh,cfg.Ns,cfg.Nd, ...
            mix1{c},mix2{c},cfg.RelTol,cfg.AbsTol);
    end
end

%% ---------------------- Chunked Monte Carlo -----------------------------
fprintf('\nRunning Monte Carlo with %.0e realizations per case...\n',cfg.Nmc);

for c = 1:nCases
    processed = 0;
    fprintf('  Case %d of %d\n',c,nCases);

    while processed < cfg.Nmc
        thisChunk = min(cfg.chunkSize,cfg.Nmc-processed);

        X = zeros(thisChunk,1);
        Y = zeros(thisChunk,1);

        if c == 4
            % Independent classical special-case validation:
            % S-R power gain is Rayleigh -> exponential with mean Omega1.
            % R-D power gain is Nakagami-m (m=2) -> Gamma(m,Omega2/m).
            mNakagami = 2;
            for n = 1:cfg.Ns
                rayleighPower = -cfg.Omega1.*log(max(rand(thisChunk,1),realmin));
                X = max(X,rayleighPower);
            end
            for n = 1:cfg.Nd
                nakagamiPower = (cfg.Omega2/mNakagami).* ...
                    randg(mNakagami.*ones(thisChunk,1));
                Y = max(Y,nakagamiPower);
            end
        else
            % General alpha-eta-mu Monte Carlo generator.
            for n = 1:cfg.Ns
                X = max(X,sampleAEM(mix1{c},thisChunk));
            end
            for n = 1:cfg.Nd
                Y = max(Y,sampleAEM(mix2{c},thisChunk));
            end
        end

        for q = 1:nSNR
            gamma1 = aVec(q).*X;
            gamma2 = bVec(q).*X.*Y;
            gammaAF = (gamma1.*gamma2)./(gamma1+gamma2+1);
            outageCounts(q,c) = outageCounts(q,c) + sum(gammaAF < gammaTh);
        end

        processed = processed + thisChunk;
        if mod(processed,1e7)==0 || processed==cfg.Nmc
            fprintf('    Processed %.0f / %.0f\n',processed,cfg.Nmc);
        end
    end
end

PoutMC = outageCounts./cfg.Nmc;

%% ----------------------- Numerical diagnostics --------------------------
fprintf('\nAnalytical / Monte Carlo comparison:\n');
for c = 1:nCases
    fprintf('\nCase %d\n',c);
    for q = 1:nSNR
        fprintf('  SNR=%2.0f dB | Exact=%.4e | MC=%.4e | events=%d\n', ...
            cfg.snr_dB(q),PoutAna(q,c),PoutMC(q,c),outageCounts(q,c));
    end
end

% Compact comparison at Pout=1e-3.
targetPout = 1e-3;
fprintf('\nRequired SNR at Pout=1e-3 (analytical interpolation):\n');
for c = 1:nCases
    snrReq = snrAtTarget(cfg.snr_dB,PoutAna(:,c),targetPout);
    fprintf('  Case %d: %.3f dB\n',c,snrReq);
end

%% ------------------------------ Plot -----------------------------------
fig = figure('Color','w','Units','inches','Position',[1 1 3.5 3.10], ...
             'PaperPositionMode','auto');
ax = axes(fig);
set(ax,'YScale','log');
hold(ax,'on'); box(ax,'on'); grid(ax,'on');

hLine = gobjects(nCases,1);
markerIndex = 1:2:nSNR;
if markerIndex(end) ~= nSNR
    markerIndex = [markerIndex,nSNR];
end

for c = 1:nCases
    hLine(c) = plot(ax,cfg.snr_dB,PoutAna(:,c),lineStyles{c}, ...
        'Color',lineColors(c,:),'LineWidth',1.1,'DisplayName',caseLabels{c});

    mcToPlot = PoutMC(markerIndex,c);
    countToPlot = outageCounts(markerIndex,c);
    mcToPlot(countToPlot < cfg.minEventsToPlot) = NaN;

    plot(ax,cfg.snr_dB(markerIndex),mcToPlot,'ko', ...
        'LineStyle','none','MarkerFaceColor','none','MarkerSize',4.5, ...
        'LineWidth',0.75,'HandleVisibility','off');
end

hSim = plot(ax,NaN,NaN,'ko','LineStyle','none','MarkerFaceColor','none', ...
    'MarkerSize',4.5,'LineWidth',0.75,'DisplayName','Simulation');

ax.YMinorGrid = 'on';
ax.XMinorGrid = 'off';
ax.FontName = 'Times New Roman';
ax.FontSize = 8;
ax.LineWidth = 0.8;
ax.TickDir = 'in';
ax.Layer = 'top';
ax.GridAlpha = 0.22;
ax.MinorGridAlpha = 0.12;

xlabel(ax,'Average source transmit SNR (dB)','Interpreter','latex','FontSize',9);
ylabel(ax,'Outage probability','Interpreter','latex','FontSize',9);
xlim(ax,[min(cfg.snr_dB) max(cfg.snr_dB)]);
ylim(ax,[1e-6 1]);

lgd = legend(ax,[hLine;hSim],[caseLabels,{'Simulation'}], ...
    'Location','southwest','Interpreter','latex');
lgd.FontName = 'Times New Roman';
lgd.FontSize = 7.2;
lgd.Box = 'on';
lgd.LineWidth = 0.8;
lgd.ItemTokenSize = [10 10];

%% ---------------------------- Save outputs ------------------------------
exportgraphics(fig,'Fig06_Fading_Generality.pdf','ContentType','vector','BackgroundColor','white');
exportgraphics(fig,'Fig06_Fading_Generality.png','Resolution',600,'BackgroundColor','white');

T = table(cfg.snr_dB(:), ...
    PoutAna(:,1),PoutMC(:,1),outageCounts(:,1), ...
    PoutAna(:,2),PoutMC(:,2),outageCounts(:,2), ...
    PoutAna(:,3),PoutMC(:,3),outageCounts(:,3), ...
    PoutAna(:,4),PoutMC(:,4),outageCounts(:,4), ...
    'VariableNames',{'SNR_dB', ...
    'Default_Exact','Default_MC','Default_Events', ...
    'Alpha15_Exact','Alpha15_MC','Alpha15_Events', ...
    'Eta02_Exact','Eta02_MC','Eta02_Events', ...
    'Rayleigh_Nakagami_Exact','Rayleigh_Nakagami_MC','Rayleigh_Nakagami_Events'});
writetable(T,'Fig06_Fading_Generality_Data.csv');

fprintf('\nSaved Fig06_Fading_Generality.pdf, Fig06_Fading_Generality.png, and Fig06_Fading_Generality_Data.csv\n');

%% =======================================================================
%                           LOCAL FUNCTIONS
% ========================================================================
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
        Fpos = zeros(size(u));

        for k = 1:numel(mix.w)
            Fpos = Fpos + mix.w(k).* ...
                   gammainc(u,mix.shape(k),'lower');
        end

        F(positive) = Fpos;
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
        zp = z(positive);
        u = mix.c.*(zp./mix.Omega).^(mix.alpha/2);

        duDz = mix.c.*(mix.alpha/2)./mix.Omega ...
             .*(zp./mix.Omega).^(mix.alpha/2-1);

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

function z = sampleAEM(mix,N)
%SAMPLEAEM Generate N alpha-eta-mu parent-gain samples.

    uCat = rand(N,1);
    cumulativeWeights = cumsum(mix.w);
    cumulativeWeights(end) = 1;

    idx = ones(N,1,'uint16');

    for k = 1:numel(cumulativeWeights)-1
        idx(uCat>cumulativeWeights(k)) = k+1;
    end

    shapePerSample = mix.shape(double(idx));
    gammaSample = randg(shapePerSample);

    z = mix.Omega.*(gammaSample./mix.c).^(2/mix.alpha);
end

function Pout = outageExactNumerical(a,b,gammaTh,Ns,Nd, ...
                                     mix1,mix2,RelTol,AbsTol)
%OUTAGEEXACTNUMERICAL Exact numerical evaluation of the one-dimensional outage representation.

    x0 = gammaTh/a;

    F1x0 = aemCDF(x0,mix1);
    L1 = F1x0^Ns;

    integrand = @(t) transformedIntegrand( ...
        t,x0,a,b,gammaTh,Ns,Nd,mix1,mix2);

    L2 = integral(integrand,0,1, ...
                  'RelTol',RelTol, ...
                  'AbsTol',AbsTol);

    Pout = min(max(L1+L2,0),1);
end

function val = transformedIntegrand(t,x0,a,b,gammaTh, ...
                                    Ns,Nd,mix1,mix2)
%TRANSFORMEDINTEGRAND Integrand after x=x0+t/(1-t).

    val = zeros(size(t));
    valid = t>0 & t<1;

    if ~any(valid)
        return;
    end

    tv = t(valid);
    x = x0 + tv./(1-tv);
    jacobian = 1./(1-tv).^2;

    denominator = b.*x.*(a.*x-gammaTh);
    yThreshold = gammaTh.*(a.*x+1)./denominator;

    F2 = aemCDF(yThreshold,mix2);
    F1 = aemCDF(x,mix1);
    f1 = aemPDF(x,mix1);

    fX = Ns.*f1.*F1.^(Ns-1);
    FY = F2.^Nd;

    temp = FY.*fX.*jacobian;
    temp(~isfinite(temp)) = 0;

    val(valid) = temp;
end

function snrValue = snrAtTarget(snr_dB,pout,target)
%SNRATTARGET Robust interpolation of the SNR at a target outage.
% Repeated outage values (for example several points equal to 1) are
% removed before calling INTERP1.

    pout = pout(:);
    snr_dB = snr_dB(:);

    valid = isfinite(pout) & pout>0 & isfinite(snr_dB);
    pout = pout(valid);
    snr_dB = snr_dB(valid);

    if isempty(pout)
        snrValue = NaN;
        return;
    end

    logP = log10(pout);

    % Sort by outage probability.
    [logPSorted,index] = sort(logP,'ascend');
    snrSorted = snr_dB(index);

    % INTERP1 requires unique sample points.
    [logPUnique,uniqueIndex] = unique(logPSorted,'stable');
    snrUnique = snrSorted(uniqueIndex);

    if numel(logPUnique)<2 || ...
       log10(target)<min(logPUnique) || ...
       log10(target)>max(logPUnique)

        snrValue = NaN;
        warning('Target outage %.1e is outside this curve range.',target);
        return;
    end

    snrValue = interp1(logPUnique,snrUnique,log10(target),'linear');
end
