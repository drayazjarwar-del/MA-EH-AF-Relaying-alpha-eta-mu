%% FIGURE 5: Relay-power saturation sensitivity
% MA-assisted EH variable-gain AF relaying over alpha-eta-mu fading.
%
% Purpose:
%   Quantify the effect of a finite relay transmit-power cap, which prevents
%   the linear-EH baseline P_R = kappa P_S X from producing unbounded relay
%   power under strong first-hop realizations.
%
% Baseline (linear EH):
%   P_R = kappa P_S X.
%
% Saturated sensitivity model:
%   P_R^sat = min(kappa P_S X, P_R,max).
%
% The script evaluates the capped model in two independent ways:
%   (i) exact one-dimensional numerical integration using the manuscript CDFs;
%   (ii) Monte Carlo simulation using the same exact Gamma-mixture sampler.
%
% Default protocol: PSR, lambda = 0.6, zeta = 0.8, R = 2 bit/s/Hz.
%
% Generated files:
%   Fig05_Relay_Power_Saturation.pdf
%   Fig05_Relay_Power_Saturation.png
%   Fig05_Relay_Power_Saturation.csv

clear; clc; close all;
rng(2026,'twister');

%% Configuration
cfg.alpha1 = 2; cfg.eta1 = 0.5; cfg.mu1 = 2; cfg.Omega1 = 1; cfg.format1 = 1;
cfg.alpha2 = 2; cfg.eta2 = 0.5; cfg.mu2 = 2; cfg.Omega2 = 1; cfg.format2 = 1;

cfg.Ns = 3;
cfg.Nd = 3;
cfg.zeta = 0.8;
cfg.lambda = 0.6;
cfg.R = 2;

cfg.snr_dB = 0:2:40;
% Normalised relay-power caps: gamma_R,max = P_R,max/N0.
cfg.cap_dB = [12 15 18];

% Increase for the final manuscript figure if desired.
cfg.Nmc = 2e7;

cfg.seriesTol = 1e-13;
cfg.Kmax = 500;
cfg.RelTol = 1e-8;
cfg.AbsTol = 1e-11;

%% Parent alpha-eta-mu distributions
mix1 = buildAEMMixture(cfg.alpha1,cfg.eta1,cfg.mu1,cfg.Omega1, ...
                       cfg.format1,cfg.seriesTol,cfg.Kmax);
mix2 = buildAEMMixture(cfg.alpha2,cfg.eta2,cfg.mu2,cfg.Omega2, ...
                       cfg.format2,cfg.seriesTol,cfg.Kmax);

fprintf('Mixture terms retained: hop 1 = %d, hop 2 = %d\n', ...
        numel(mix1.w),numel(mix2.w));

%% Generate selected gains once and reuse across all SNR/cap cases
fprintf('Generating %g selected first-hop gains...\n',cfg.Nmc);
X = zeros(cfg.Nmc,1);
for n = 1:cfg.Ns
    X = max(X,sampleAEM(mix1,cfg.Nmc));
end

fprintf('Generating %g selected second-hop gains...\n',cfg.Nmc);
Y = zeros(cfg.Nmc,1);
for m = 1:cfg.Nd
    Y = max(Y,sampleAEM(mix2,cfg.Nmc));
end

%% PSR parameters
snrLin = 10.^(cfg.snr_dB/10);
gammaTh = 2^(2*cfg.R)-1;

nSNR = numel(snrLin);
nCap = numel(cfg.cap_dB);

PoutLinearAna = zeros(nSNR,1);
PoutLinearMC  = zeros(nSNR,1);
PoutCapAna = zeros(nSNR,nCap);
PoutCapMC  = zeros(nSNR,nCap);

for q = 1:nSNR
    a = (1-cfg.lambda)*snrLin(q);
    b = cfg.zeta*cfg.lambda*snrLin(q);

    % Linear-EH analytical baseline.
    PoutLinearAna(q) = outageExactNumericalLinear( ...
        a,b,gammaTh,cfg.Ns,cfg.Nd,mix1,mix2,cfg.RelTol,cfg.AbsTol);

    gamma1 = a.*X;
    gamma2 = b.*X.*Y;
    gammaAF = (gamma1.*gamma2)./(gamma1+gamma2+1);
    PoutLinearMC(q) = mean(gammaAF < gammaTh);

    for cidx = 1:nCap
        gammaRmax = 10^(cfg.cap_dB(cidx)/10); % P_R,max/N0

        % Exact numerical evaluation under the capped relay-power model.
        PoutCapAna(q,cidx) = outageExactNumericalCapped( ...
            a,b,gammaRmax,gammaTh,cfg.Ns,cfg.Nd, ...
            mix1,mix2,cfg.RelTol,cfg.AbsTol);

        % Monte Carlo evaluation under the same capped model.
        gamma2sat = min(b.*X,gammaRmax).*Y;
        gammaAFsat = (gamma1.*gamma2sat)./(gamma1+gamma2sat+1);
        PoutCapMC(q,cidx) = mean(gammaAFsat < gammaTh);
    end
end

%% High-SNR floors under relay-power saturation
fprintf('\nHigh-SNR floors for the capped model:\n');
for cidx = 1:nCap
    gammaRmax = 10^(cfg.cap_dB(cidx)/10);
    floorVal = aemCDF(gammaTh/gammaRmax,mix2)^cfg.Nd;
    fprintf('  cap = %2.0f dB : floor = %.6e\n',cfg.cap_dB(cidx),floorVal);
end

%% Plot
fig = figure('Color','w','Units','inches','Position',[1 1 3.5 3.10]);
ax = axes(fig); hold(ax,'on'); set(ax,'YScale','log');

% Distinct line styles; markers are reserved for Monte Carlo.
h0 = plot(ax,cfg.snr_dB,PoutLinearAna,'-','LineWidth',1.2, ...
          'DisplayName','Linear EH (exact)');

styles = {'--','-.',':'};
hCap = gobjects(nCap,1);
for cidx = 1:nCap
    hCap(cidx) = plot(ax,cfg.snr_dB,PoutCapAna(:,cidx),styles{cidx}, ...
        'LineWidth',1.2, ...
        'DisplayName',sprintf('$\\bar{\\gamma}_{R,\\max}=%g$ dB',cfg.cap_dB(cidx)));
end

% Simulation circles. Plot every other SNR point to avoid clutter.
markIdx = 1:2:nSNR;
plot(ax,cfg.snr_dB(markIdx),replaceZero(PoutLinearMC(markIdx)),'o', ...
     'LineStyle','none','Color','k','MarkerFaceColor','none', ...
     'MarkerSize',4.5,'LineWidth',0.8,'HandleVisibility','off');
for cidx = 1:nCap
    plot(ax,cfg.snr_dB(markIdx),replaceZero(PoutCapMC(markIdx,cidx)),'o', ...
         'LineStyle','none','Color','k','MarkerFaceColor','none', ...
         'MarkerSize',4.5,'LineWidth',0.8,'HandleVisibility','off');
end
hSim = plot(ax,NaN,NaN,'ko','LineStyle','none','MarkerFaceColor','none', ...
            'MarkerSize',4.5,'LineWidth',0.8,'DisplayName','Simulation');

grid(ax,'on'); box(ax,'on');
ax.YMinorGrid = 'on'; ax.XMinorGrid = 'off';
ax.GridAlpha = 0.22; ax.MinorGridAlpha = 0.12;
ax.FontName = 'Times New Roman'; ax.FontSize = 8;
ax.LineWidth = 0.8; ax.TickDir = 'in'; ax.Layer = 'top';

xlabel(ax,'Average source transmit SNR (dB)','FontSize',9);
ylabel(ax,'Outage probability','FontSize',9);

legend(ax,[h0;hCap;hSim],'Location','southwest','Interpreter','latex', ...
       'FontName','Times New Roman','FontSize',8,'Box','on');

xlim(ax,[min(cfg.snr_dB) max(cfg.snr_dB)]);
ylim(ax,[1e-7 1]);

exportgraphics(fig,'Fig05_Relay_Power_Saturation.pdf', ...
               'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,'Fig05_Relay_Power_Saturation.png', ...
               'Resolution',600);

%% Save data
T = table(cfg.snr_dB(:),PoutLinearAna,PoutLinearMC, ...
          PoutCapAna(:,1),PoutCapMC(:,1), ...
          PoutCapAna(:,2),PoutCapMC(:,2), ...
          PoutCapAna(:,3),PoutCapMC(:,3), ...
    'VariableNames',{'SNR_dB','Linear_Exact','Linear_MC', ...
      'Cap12dB_Exact','Cap12dB_MC','Cap15dB_Exact','Cap15dB_MC', ...
      'Cap18dB_Exact','Cap18dB_MC'});
writetable(T,'Fig05_Relay_Power_Saturation.csv');

fprintf('\nFigure and data exported.\n');

%% ======================================================================
% Local functions
% =======================================================================
function mix = buildAEMMixture(alpha,eta,mu,Omega,formatNo,tol,Kmax)
    [h,H] = etaAuxiliary(eta,formatNo);
    c = 2*mu*h;
    if abs(H) < 1e-15
        k = 0; w = 1;
    else
        k = (0:Kmax).';
        logw = log(2)+0.5*log(pi)+2*mu*log(mu)+mu*log(h)-gammaln(mu) ...
             +2*k.*log(abs(mu*H))-gammaln(k+1)-gammaln(mu+k+0.5) ...
             -2*(mu+k).*log(c)+gammaln(2*(mu+k));
        wAll = exp(logw);
        cumulativeMass = cumsum(wAll);
        last = find(cumulativeMass >= 1-tol,1,'first');
        if isempty(last), last = numel(wAll); end
        k = k(1:last); w = wAll(1:last); w = w/sum(w);
    end
    mix.alpha=alpha; mix.eta=eta; mix.mu=mu; mix.Omega=Omega;
    mix.formatNo=formatNo; mix.h=h; mix.H=H; mix.c=c;
    mix.k=k; mix.shape=2*(mu+k); mix.w=w;
end

function [h,H] = etaAuxiliary(eta,formatNo)
    switch formatNo
        case 1
            h=(1+eta)^2/(4*eta); H=(1-eta^2)/(4*eta);
        case 2
            h=1/(1-eta^2); H=eta/(1-eta^2);
        otherwise
            error('formatNo must be 1 or 2.');
    end
end

function F = aemCDF(z,mix)
    originalSize=size(z); z=z(:).'; F=zeros(size(z));
    positive=z>0 & isfinite(z);
    if any(positive)
        u=mix.c.*(z(positive)./mix.Omega).^(mix.alpha/2);
        Fp=zeros(size(u));
        for kk=1:numel(mix.w)
            Fp=Fp+mix.w(kk).*gammainc(u,mix.shape(kk),'lower');
        end
        F(positive)=Fp;
    end
    F(isinf(z)&z>0)=1; F(z<=0)=0; F=min(max(F,0),1);
    F=reshape(F,originalSize);
end

function f = aemPDF(z,mix)
    originalSize=size(z); z=z(:).'; f=zeros(size(z));
    positive=z>0 & isfinite(z);
    if any(positive)
        zp=z(positive);
        u=mix.c.*(zp./mix.Omega).^(mix.alpha/2);
        duDz=mix.c.*(mix.alpha/2)./mix.Omega.*(zp./mix.Omega).^(mix.alpha/2-1);
        fp=zeros(size(zp));
        for kk=1:numel(mix.w)
            s=mix.shape(kk);
            gammaPDF=exp((s-1).*log(u)-u-gammaln(s));
            fp=fp+mix.w(kk).*gammaPDF.*duDz;
        end
        f(positive)=fp;
    end
    f(z<=0)=0; f=reshape(f,originalSize);
end

function z = sampleAEM(mix,N)
    uCat=rand(N,1); cw=cumsum(mix.w); cw(end)=1;
    idx=ones(N,1,'uint16');
    for kk=1:numel(cw)-1
        idx(uCat>cw(kk))=kk+1;
    end
    shapePerSample=mix.shape(double(idx));
    gammaSample=randg(shapePerSample);
    z=mix.Omega.*(gammaSample./mix.c).^(2/mix.alpha);
end

function Pout = outageExactNumericalLinear(a,b,gammaTh,Ns,Nd,mix1,mix2,RelTol,AbsTol)
    x0=gammaTh/a;
    L1=aemCDF(x0,mix1)^Ns;
    integrand=@(t) transformedIntegrandLinear(t,x0,a,b,gammaTh,Ns,Nd,mix1,mix2);
    L2=integral(integrand,0,1,'RelTol',RelTol,'AbsTol',AbsTol);
    Pout=min(max(L1+L2,0),1);
end

function val = transformedIntegrandLinear(t,x0,a,b,gammaTh,Ns,Nd,mix1,mix2)
    val=zeros(size(t)); valid=t>0 & t<1; if ~any(valid), return; end
    tv=t(valid); x=x0+tv./(1-tv); jac=1./(1-tv).^2;
    yth=gammaTh.*(a.*x+1)./(b.*x.*(a.*x-gammaTh));
    F2=aemCDF(yth,mix2); F1=aemCDF(x,mix1); f1=aemPDF(x,mix1);
    fX=Ns.*f1.*F1.^(Ns-1); FY=F2.^Nd;
    tmp=FY.*fX.*jac; tmp(~isfinite(tmp))=0; val(valid)=tmp;
end

function Pout = outageExactNumericalCapped(a,b,gammaRmax,gammaTh,Ns,Nd,mix1,mix2,RelTol,AbsTol)
    x0=gammaTh/a;
    L1=aemCDF(x0,mix1)^Ns;
    integrand=@(t) transformedIntegrandCapped(t,x0,a,b,gammaRmax,gammaTh,Ns,Nd,mix1,mix2);
    L2=integral(integrand,0,1,'RelTol',RelTol,'AbsTol',AbsTol);
    Pout=min(max(L1+L2,0),1);
end

function val = transformedIntegrandCapped(t,x0,a,b,gammaRmax,gammaTh,Ns,Nd,mix1,mix2)
    val=zeros(size(t)); valid=t>0 & t<1; if ~any(valid), return; end
    tv=t(valid); x=x0+tv./(1-tv); jac=1./(1-tv).^2;
    relaySNR=min(b.*x,gammaRmax); % P_R^sat/N0
    yth=gammaTh.*(a.*x+1)./(relaySNR.*(a.*x-gammaTh));
    F2=aemCDF(yth,mix2); F1=aemCDF(x,mix1); f1=aemPDF(x,mix1);
    fX=Ns.*f1.*F1.^(Ns-1); FY=F2.^Nd;
    tmp=FY.*fX.*jac; tmp(~isfinite(tmp))=0; val(valid)=tmp;
end

function y=replaceZero(y)
    y(y<=0)=NaN;
end
