%% Generate Tumor Drug profile for BRAFi/MEKi pairs in a avascular spherical tumor
%% Miles Miller, Thomas Ng

function [tmat, tsum, tfrac, nDB, NP, BSITES_T] = simulate_TUMOR(input);
%% Parameter values 
% All volumes in mL
% All time in min
% DOSE in mg
% Lengths in cm (cm^3 = mL)

close all

global NSpecies;
global DOSESTD;
global nDB NP;
global BSITES BSITES_T;
global eps_d eps_a;

NPARAM = 30;
BSITES_T = 1e-9;                                            % molarity of target protein % PMC5293153 [mol ml-1]

%% Input parameters
RCAP            = input.RCAP;                               %  capillary radius 15um
q               = input.q;                                  % absorption rate
deltaX          = input.deltaX;                             % integrating distance within the tumor in cm 

ATP_conc        = input.ATP_conc;                           % Intracellular ATP concentration
k_unbind_ATP    = input.k_unbind_ATP;                       % ATP-BRAF dissociation
k_bind_ATP      = input.k_bind_ATP;                         % ATP-BRAF association
k_cat_ATP       = input.k_cat_ATP;                          % catabolism of ATP-bound moieties

DOSEmgkg        = input.DOSEmgkg;                           % Dose administered in mg
BIO             = input.BIO;                                % Oral bioavailability (fraction)
Mass            = input.Mass;                               % Mass of subject (Set as 1 for human, mouse mass in kg)

V1      = input.V1;                                        % volume of distribution (mL)  

MW              = input.MW;                                 % MW in g/mol


TimeLen         = input.TimeLen;                            % Length of simulation (min)
DoseInt         = input.DoseInt;                            % Dosing interval (min)

k_bind          = input.k_bind;                             % drug target association
k_unbind        = input.k_unbind;                           % drug target dissociation

k_el1           = input.k_el1;                              % fast phase plasma elimination (1/min)
k_el2           = input.k_el2;                              % slow phase plasma elimination (1/min))
k_elA           = input.k_elA;                              % slow phase contribution (fraction)

kPerm1          = input.kPerm1;                             % permeability of tumor vasculature to unbound drug
Deff            = input.Deff;                               % effective diffusion of unbound drug


k_uptake  = input.k_uptake;                                 % partitioning into neutral lipid min-1
k_release = input.k_release;                                % partitioning into aqueous (lipid release)

DeffA     = input.DeffA;                                    % Effective diffusion of albumin bound drug
kPerm1A   = input.kPerm1A;                                  % permeability of tumor vasculature to albumin
konALB    = input.konALB;                                   % albumin-drug association
koffALB   = input.koffALB;                                  % albumin drug dissociation, calculated from KD and known drug plasma protein binding

k_uptakeALB = input.k_uptakeALB;                            % cellular drug-albumin uptake
FF = input.FF;                                              % fraction of neutral lipid
k_catabALB = input.k_catabALB;                              % drug release from cellular albumin-drug complex

nDB = input.nDB;                                            % Rmax max radius of model, 500um

%% DO NOT MODIFY BELOW

KpT       = 1; %16  blood tissue partition coefficient
Kp_norm   = 0;%18  300; %partition coefficient whole body
PermAnorm = 0;%20  5e-3; %%an effective P*S term...need to update
Kp_ALBnorm = 1; %%21  albumin is >3x concentrated in serum than tissue
SAtumor = 2 * pi * RCAP * deltaX;  %6  surface area of tumor vessel
Vt  = pi*(RCAP+deltaX).^2*deltaX - pi*(RCAP)^2*deltaX; %7  peritumor volume
Vs  = 0; %8   
DOSEmgkg = DOSEmgkg*BIO;
MOUSEmass = Mass;%0.02; % 11  pin(21) kg;
DOSE      = DOSEmgkg .* MOUSEmass .* 1e-3 ./ MW; % 12  dose in moles

NP = 4; %species before the method of lines section (unbound blood, stomach, bound blood)
NSpecies = NP + nDB*8;
%%

pin = [q V1 k_el1 ...
    kPerm1 SAtumor ...
    Vt Vs RCAP ...
    DOSEmgkg MOUSEmass DOSE ...
    Deff deltaX ...
    k_uptake k_release KpT ...
    k_el2 Kp_norm k_elA ...
    PermAnorm Kp_ALBnorm DeffA ...
    kPerm1A konALB koffALB k_uptakeALB FF ...
    k_catabALB k_bind k_unbind ...
    k_bind_ATP k_unbind_ATP, ATP_conc, k_cat_ATP, TimeLen, DoseInt];

x0 = log10(pin);

[cost, tmat,tsum] = main_solveODE(x0);
clear tplot tfrac;

for plotIndSamples = 1 : 0
    
    titleind = {'Free'
        'lipid binding sites'
        'lipid bound'
        'alb bound'
        'cell alb bound'
        'free target'
        'drug bound target'};
    for iii = 0 : 6
        for nn = 1 : size(tmat,1)
            tplot(nn,:)=tmat(nn,NP+nDB*iii+1:NP+nDB*(iii+1));
        end
        
        RNN = 1:size(tplot,1);
        CNN = 1:size(tplot,2);
        NMt = length(RNN);
        NCt = length(CNN);
        ztt = zeros(NMt,1);
        pmats = [tplot ztt];
        ztt = zeros(1,NCt+1);
        pmats = [pmats; ztt];
        
        figure; hold all;
        h=pcolor(log10(pmats));
        colormap(jet);
        set(h,'edgecolor','none');
        %     set(gca,'clim',[-10 -5.7]);
        axis xy;
        ylim([1 NMt+1]);
        xlim([1 NCt+1]);
        set(gca,'xtick',[1.5:6:NCt+.5]);
        xtickangle(45);
        set(gca,'ytick',[1.5:10:NMt+.5]);
        title(titleind{iii+1});
        colorbar;
    end
    
end

for plotCost = 1
    
    RNN = 1:size(tsum,1);
    CNN = 1:size(tsum,2);
    NMt = length(RNN);
    NCt = length(CNN);
    ztt = zeros(NMt,1);
    pmats = [tsum ztt];
    ztt = zeros(1,NCt+1);
    pmats = [pmats; ztt];
    
    figure; hold all;
    h=pcolor((pmats));
    colormap(jet);
    set(h,'edgecolor','none');
%     set(gca,'clim',[-10 -5.7]);
    axis xy;
    ylim([1 NMt+1]);
    xlim([1 NCt+1]);
    set(gca,'xtick',[1.5:6:NCt+.5]);
    % set(gca,'xticklabel',xlab(CNN));
    xtickangle(45);
    set(gca,'ytick',[1.5:10:NMt+.5]);
    % set(gca,'yticklabel',DVEC);
    title('sum drug');
    xlabel('Distance from tumor core')
    ylabel('Time')
    colorbar;
   
    figure;
    contour((tsum));

    
    for nn = 1 : size(tmat,1)
        tfrac(nn,:)=tmat(nn,NP+nDB*6+1:NP+nDB*(6+1))./BSITES_T./1e3;
    end
    
    RNN = 1:size(tfrac,1);
    CNN = 1:size(tfrac,2);
    NMt = length(RNN);
    NCt = length(CNN);
    ztt = zeros(NMt,1);
    pmats = [tfrac ztt];
    ztt = zeros(1,NCt+1);
    pmats = [pmats; ztt];
    
    figure; hold all;
    h=pcolor((pmats));
    colormap(jet);
    set(h,'edgecolor','none');
    set(gca,'clim',[0 1]);
    axis xy;
    ylim([1 NMt+1]);
    xlim([1 NCt+1]);
    set(gca,'xtick',[1.5:6:NCt+.5]);
    xtickangle(45);
    set(gca,'ytick',[1.5:10:NMt+.5]);
    title('frac occ');
    xlabel('Distance from tumor core')
    ylabel('Time')
    colorbar;
    
    figure;
    contour(tfrac,[0.1:0.1:1],'linewidth',1);
   
end

for plotTestFractions = 1:0
    
    for nn = 1 : size(tmat,1)
        tfrac(nn,:)=tmat(nn,NP+1:NP+nDB)./tmat(nn,NP+nDB*2+1:NP+nDB*3);
    end
    
    RNN = 1:size(tplot,1);
    CNN = 1:size(tplot,2);
    NMt = length(RNN);
    NCt = length(CNN);
    ztt = zeros(NMt,1);
    pmats = [tfrac ztt];
    ztt = zeros(1,NCt+1);
    pmats = [pmats; ztt];
    
    figure; hold all;
    h=pcolor(log10(pmats));
    colormap(jet);
    set(h,'edgecolor','none');
    %     set(gca,'clim',[-10 -5.7]);
    axis xy;
    ylim([1 NMt+1]);
    xlim([1 NCt+1]);
    set(gca,'xtick',[1.5:6:NCt+.5]);
    xtickangle(45);
    set(gca,'ytick',[1.5:10:NMt+.5]);
    title('frac lipid');
    colorbar;
    
    
    for nn = 1 : size(tmat,1)
        tfrac(nn,:)=tmat(nn,NP+1:NP+nDB)./tmat(nn,NP+nDB*3+1:NP+nDB*4);
    end
    
    RNN = 1:size(tplot,1);
    CNN = 1:size(tplot,2);
    NMt = length(RNN);
    NCt = length(CNN);
    ztt = zeros(NMt,1);
    pmats = [tfrac ztt];
    ztt = zeros(1,NCt+1);
    pmats = [pmats; ztt];
    
    figure; hold all;
    h=pcolor(log10(tfrac));
    colormap(jet);
    set(h,'edgecolor','none');
    %     set(gca,'clim',[-10 -5.7]);
    axis xy;
    ylim([1 NMt+1]);
    xlim([1 NCt+1]);
    set(gca,'xtick',[1.5:6:NCt+.5]);
    xtickangle(45);
    set(gca,'ytick',[1.5:10:NMt+.5]);
    title('frac alb bound');
    colorbar;
    
end

for plotLines = 1
    
   Tvec2 = [1/60 1 4 23]*60;
%     Tvec2 = [1/60 2 23 26 47 50 71 74 119 122 239 242 479 482]*60;
  %  Tvec2 = [1/60 2 23 26 47 50 71 74 119 122]*60;
    cmat = jet(length(Tvec2));
    figure; hold all;
    for nn = 1 : length(Tvec2)
        plot(1:nDB,tsum(Tvec2(nn),:),'Color',cmat(nn,:),'linewidth',2);
    end
    
    title('Total drug');
     legend('5 min', '1 hr','4 hr', '23 hr');
%    legend('5 min', '1 hr', '1 day', '1day+1hr','2 day','2 day+1hr','3 day','3 day+1hr','5 day','5 day+1hr');
    
end

for COMPCOST = 1
    clear costy;
    costy(1,:) = [max(tfrac(:)) mean(tfrac(:)) std(tfrac(:))];
    costy(2,:) = [max(tfrac(:,1)) mean(tfrac(:,1)) std(tfrac(:,1))];
    costy(3,:) = [max(tfrac(:,length(tfrac(1,:)))) mean(tfrac(:,length(tfrac(1,:)))) std(tfrac(:,length(tfrac(1,:))))];

    costy(4,:) = [max(tsum(:)) mean(tsum(:)) std(tsum(:))];
    costy(5,:) = [max(tsum(:,1)) mean(tsum(:,1)) std(tsum(:,1))];
    costy(6,:) = [max(tsum(:,length(tsum(1,:)))) mean(tsum(:,length(tsum(1,:)))) std(tsum(:,length(tsum(1,:))))];
    
    
end

  for nn = 1 : size(tmat,1)
        tfrac(nn,:)=tmat(nn,NP+nDB*6+1:NP+nDB*(6+1))./BSITES_T./1e3;
  end
