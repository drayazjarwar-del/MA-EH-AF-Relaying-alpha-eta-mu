%% FIGURE 3: Effect of movable-antenna deployment under PSR
% Movable-antenna-assisted energy-harvesting AF relaying over
% alpha-eta-mu fading channels.
%
% Compared configurations:
%   1) Ns = 1, Nd = 1 : fixed antennas at source and destination
%   2) Ns = 3, Nd = 1 : movable antenna only at the source
%   3) Ns = 1, Nd = 3 : movable antenna only at the destination
%   4) Ns = 3, Nd = 3 : movable antennas at both terminals
%
% Analytical curves:
%   Exact numerical evaluation of the one-dimensional outage representation.
%
% Simulation markers:
%   Monte Carlo evaluation of the end-to-end AF SNR.
%
% The code uses chunked Monte Carlo simulation so that 1e7 realizations
% can be evaluated without storing all channel samples in memory.
%
% MATLAB requirement:
%   Local functions in scripts require MATLAB R2016b or newer.
%
% Generated files:
%   Fig03_MA_Deployment_PSR.pdf   (vector PDF for LaTeX)
%   Fig03_MA_Deployment_PSR.png   (600-dpi backup)
%   Fig03_MA_Deployment_PSR.eps
%   Fig03_MA_Deployment_PSR_Data.csv

clear; clc; close all;
rng(2026,'twister');

set(groot,'defaultFigureVisible','on');

% Open a visible status window immediately so the user can confirm
% that the script has started.
statusFig = figure('Color','w', ...
                   'Visible','on', ...
                   'Name','Figure 3 calculation status', ...
                   'NumberTitle','off', ...
                   'Position',[300 300 560 190]);

statusAx = axes(statusFig);
axis(statusAx,'off');

statusText = text(statusAx,0.5,0.60, ...
    'Figure 3 calculation is running...', ...
    'HorizontalAlignment','center', ...
    'FontSize',13);

statusSubText = text(statusAx,0.5,0.37, ...
    'Please wait while the analytical curves are evaluated.', ...
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

cfg.zeta   = 0.8;
cfg.lambda = 0.6;                   % PSR power-splitting factor
cfg.R      = 2;                     % bit/s/Hz

cfg.snr_dB = 0:2:30;                % average source transmit SNR
cfg.Nmc    = 1e8;                   % journal-quality Monte Carlo samples
cfg.chunkSize = 2e8;                % memory-efficient simulation block

cfg.seriesTol = 1e-13;
cfg.Kmax      = 500;

cfg.RelTol = 1e-7;
cfg.AbsTol = 1e-11;

% Antenna-deployment configurations: [Ns, Nd]
configs = [1 1;
           3 1;
           1 3;
           3 3];

configLabels = { ...
 '$N_s=1,\;N_d=1$', ...
    '$N_s=3,\;N_d=1$', ...
    '$N_s=1,\;N_d=3$', ...
    '$N_s=3,\;N_d=3$'};

% Analytical curves use distinct colors while preserving line styles.
lineStyles = {':','--','-.','-'};

lineColors = [
    0       0.4470 0.7410;   % Ns=1, Nd=1  blue
    0.8500  0.3250 0.0980;   % Ns=3, Nd=1  orange/red
    0.4660  0.6740 0.1880;   % Ns=1, Nd=3  green
    0.4940  0.1840 0.5560    % Ns=3, Nd=3  purple
];

%% -------------------- Build parent fading distributions -----------------
mix1 = buildAEMMixture(cfg.alpha1, cfg.eta1, cfg.mu1, ...
                       cfg.Omega1, cfg.format1, ...
                       cfg.seriesTol, cfg.Kmax);

mix2 = buildAEMMixture(cfg.alpha2, cfg.eta2, cfg.mu2, ...
                       cfg.Omega2, cfg.format2, ...
                       cfg.seriesTol, cfg.Kmax);

fprintf('First-hop mixture terms retained : %d\n', numel(mix1.w));
fprintf('Second-hop mixture terms retained: %d\n', numel(mix2.w));

%% ---------------------- Basic PDF validation ----------------------------
pdfMass1 = integral(@(z) aemPDF(z,mix1),0,Inf, ...
                    'RelTol',1e-8,'AbsTol',1e-11);
pdfMass2 = integral(@(z) aemPDF(z,mix2),0,Inf, ...
                    'RelTol',1e-8,'AbsTol',1e-11);

fprintf('Integral of first-hop parent PDF : %.12f\n',pdfMass1);
fprintf('Integral of second-hop parent PDF: %.12f\n\n',pdfMass2);

if abs(pdfMass1-1)>1e-6 || abs(pdfMass2-1)>1e-6
    warning(['A parent PDF does not integrate sufficiently close to one. ', ...
             'Increase cfg.Kmax or inspect the fading parameters.']);
end

%% ---------------------- PSR protocol parameters -------------------------
snrLin  = 10.^(cfg.snr_dB/10);
gammaTh = 2^(2*cfg.R)-1;

% Under PSR:
%   gamma1 = aX,  a = (1-lambda) Ps/N0
%   gamma2 = bXY, b = zeta lambda Ps/N0
aVec = (1-cfg.lambda).*snrLin;
bVec = cfg.zeta*cfg.lambda.*snrLin;

nSNR    = numel(snrLin);
nConfig = size(configs,1);

%% ---------------- Exact analytical outage probability ------------------
PoutAna = zeros(nSNR,nConfig);

fprintf('Evaluating exact analytical curves...\n');
statusText.String = 'Evaluating exact analytical curves...';
statusSubText.String = 'The final graph will appear after all calculations finish.';
drawnow;

for c = 1:nConfig
    Ns = configs(c,1);
    Nd = configs(c,2);

    fprintf('  Configuration Ns=%d, Nd=%d\n',Ns,Nd);
    if isgraphics(statusText)
        statusText.String = sprintf('Analytical: (N_s,N_d)=(%d,%d)',Ns,Nd);
        drawnow;
    end

    for q = 1:nSNR
        PoutAna(q,c) = outageExactNumerical( ...
            aVec(q),bVec(q),gammaTh,Ns,Nd, ...
            mix1,mix2,cfg.RelTol,cfg.AbsTol);
    end
end

fprintf('\n');

%% ---------------- Chunked Monte Carlo simulation -----------------------
outageCounts = zeros(nSNR,nConfig);
processed = 0;
maxNs = max(configs(:,1));
maxNd = max(configs(:,2));

fprintf('Running Monte Carlo simulation with %.0e realizations...\n',cfg.Nmc);
if isgraphics(statusText)
    statusText.String = 'Running Monte Carlo simulation...';
    statusSubText.String = sprintf('0 of %.0f realizations processed',cfg.Nmc);
    drawnow;
end

while processed < cfg.Nmc
    thisChunk = min(cfg.chunkSize,cfg.Nmc-processed);

    % Generate first-hop candidate gains.
    Xcand = zeros(thisChunk,maxNs);
    for n = 1:maxNs
        Xcand(:,n) = sampleAEM(mix1,thisChunk);
    end

    % Generate second-hop candidate gains.
    Ycand = zeros(thisChunk,maxNd);
    for m = 1:maxNd
        Ycand(:,m) = sampleAEM(mix2,thisChunk);
    end

    % Selected gains required by the four configurations.
    X1 = Xcand(:,1);
    X3 = max(Xcand(:,1:3),[],2);
    Y1 = Ycand(:,1);
    Y3 = max(Ycand(:,1:3),[],2);

    for c = 1:nConfig
        Ns = configs(c,1);
        Nd = configs(c,2);

        if Ns == 1
            X = X1;
        else
            X = X3;
        end

        if Nd == 1
            Y = Y1;
        else
            Y = Y3;
        end

        for q = 1:nSNR
            gamma1  = aVec(q).*X;
            gamma2  = bVec(q).*X.*Y;
            gammaAF = (gamma1.*gamma2)./(gamma1+gamma2+1);

            outageCounts(q,c) = outageCounts(q,c) + ...
                                sum(gammaAF < gammaTh);
        end
    end

    processed = processed + thisChunk;

    if mod(processed,1e6)==0 || processed==cfg.Nmc
        fprintf('  Processed %.0f of %.0f realizations\n', ...
                processed,cfg.Nmc);
    end
end

PoutMC = outageCounts./cfg.Nmc;

%% ----------------------- Display numerical agreement -------------------
fprintf('\nAnalytical and simulation results:\n');

for c = 1:nConfig
    fprintf('\nConfiguration Ns=%d, Nd=%d\n', ...
            configs(c,1),configs(c,2));

    for q = 1:nSNR
        fprintf(['  SNR=%4.1f dB | Analytical=%.4e | ', ...
                 'Simulation=%.4e | outages=%d\n'], ...
                 cfg.snr_dB(q),PoutAna(q,c),PoutMC(q,c), ...
                 outageCounts(q,c));
    end
end

%% ------------------------------ Plot -----------------------------------
if isgraphics(statusFig)
    close(statusFig);
end

fig = figure('Color','w', ...
             'Units','inches', ...
             'Position',[1 1 3.5 3.10], ...
             'Visible','on', ...
             'Renderer','opengl');



ax = axes(fig);
set(ax,'YScale','log');
hold(ax,'on');

hLine = gobjects(nConfig,1);

% Use alternate SNR points for simulation markers to reduce clutter.
markerIndex = 1:2:nSNR;
if markerIndex(end) ~= nSNR
    markerIndex = [markerIndex nSNR];
end

for c = 1:nConfig
    % Analytical lines: all black.
    hLine(c) = plot(ax,cfg.snr_dB,PoutAna(:,c), ...
        lineStyles{c}, ...
        'Color',lineColors(c,:), ...
        'LineWidth',1, ...
        'DisplayName',configLabels{c});

    % Simulation markers: common black open circles.
    mcToPlot = PoutMC(markerIndex,c);
    mcToPlot(mcToPlot==0) = NaN;

    plot(ax,cfg.snr_dB(markerIndex),mcToPlot,'o', ...
        'LineStyle','none', ...
        'Color','k', ...
        'MarkerFaceColor','none',...
        'MarkerSize',5, ...
        'LineWidth',0.8, ...
        'HandleVisibility','off');
end

% Dummy marker for the simulation legend entry.
hSim = plot(ax,NaN,NaN,'ko', ...
    'LineStyle','none', ...
    'MarkerFaceColor','w', ...
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

xlabel(ax, ...
    'Average source transmit SNR', ...
    'FontSize',9);

ylabel(ax, ...
    'Outage probability', ...
    'FontSize',9);

legendLabels = { ...
    '$N_s=1,\;N_d=1$', ...
    '$N_s=3,\;N_d=1$', ...
    '$N_s=1,\;N_d=3$', ...
    '$N_s=3,\;N_d=3$', ...
    'Simulation'};

lgd = legend(ax,[hLine;hSim], ...
    legendLabels, ...
    'Location','southwest', ...
    'Interpreter','latex');

lgd.FontSize = 11;      % <-- increase from 8 to 12
lgd.Box = 'on';
lgd.LineWidth = 1;
lgd.ItemTokenSize = [11 11];

drawnow
lgd.Position(3) = lgd.Position(3)*1.15;

xlim(ax,[min(cfg.snr_dB) max(cfg.snr_dB)]);
ylim(ax,[1e-6 1]);

xticks(ax,0:5:30);
yticks(ax,10.^(-6:0));

drawnow;
figure(fig);
shg;

%% ---------------------------- Save outputs ------------------------------
set(fig,'PaperPositionMode','auto');

print(fig,'Fig03_MA_Deployment_PSR.pdf','-dpdf','-painters');

exportgraphics(fig,'Fig03_MA_Deployment_PSR.png', ...
               'Resolution',600, ...
               'BackgroundColor','white');

print(fig,'Fig03_MA_Deployment_PSR.eps','-depsc2','-r600');

T = table(cfg.snr_dB(:), ...
          PoutAna(:,1),PoutMC(:,1),outageCounts(:,1), ...
          PoutAna(:,2),PoutMC(:,2),outageCounts(:,2), ...
          PoutAna(:,3),PoutMC(:,3),outageCounts(:,3), ...
          PoutAna(:,4),PoutMC(:,4),outageCounts(:,4), ...
    'VariableNames', { ...
        'SNR_dB', ...
        'Ns1_Nd1_Analytical','Ns1_Nd1_MC','Ns1_Nd1_OutageCount', ...
        'Ns3_Nd1_Analytical','Ns3_Nd1_MC','Ns3_Nd1_OutageCount', ...
        'Ns1_Nd3_Analytical','Ns1_Nd3_MC','Ns1_Nd3_OutageCount', ...
        'Ns3_Nd3_Analytical','Ns3_Nd3_MC','Ns3_Nd3_OutageCount'});

writetable(T,'Fig03_MA_Deployment_PSR_Data.csv');

fprintf('\nFigure 2 and numerical results have been exported.\n');

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
