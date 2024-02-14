function [cost, tmat,tsum] = main_solveODE(x0);

global NSpecies;
global DOSESTD;
global nDB NP;
global BSITES BSITES_T;
global eps_d eps_a;


%% Parameter values 
% All volumes in mL
% All time in min
% DOSE in mg

%% Parameter values 
piny = 10.^x0;


% initialize parameters (explained elsewhere)
for initializeparam = 1

q       = piny(1);
V1      = piny(2);
k_el1   = piny(3);

kPerm1  = piny(4);
SAtumor = piny(5);

Vt  = piny(6);
Vs  = piny(7);
RCAP = piny(8);

DOSEmgkg  = piny(9);
MOUSEmass = piny(10);
DOSE      = piny(11);

Deff      = piny(12);
deltaX    = piny(13);

k_uptake  = piny(14);
k_release = piny(15);

KpT       = piny(16);

k_el2 = piny(17);
Kp_norm = piny(18);
k_elA = piny(19);
PermAnorm = piny(20);
Kp_ALBnorm = piny(21);
DeffA = piny(22);
kPerm1A = piny(23);
konALB = piny(24);
koffALB = piny(25);
k_uptakeALB = piny(26);
FF = piny(27);


k_catabALB = piny(28);
k_bind     = piny(29);
k_unbind   = piny(30);
k_bind_ATP = piny(31);
k_unbind_ATP=piny(32);
ATP_conc   = piny(33);
k_cat_ATP  = piny(34);
TimeLen   = piny(35);

end



% initial simulation concentrations
y0 = zeros(NSpecies,1);
y0(NP+nDB*5+1:NP+nDB*6) = BSITES_T; %Free target binding sites % PO


%%
% epsilon - drug
eps_d = 0.4; % similar to doxorubicin as model sm molecule
eps_a = 0.2; % extracellular 

%% Solve the systems of ODEs
options = odeset('AbsTol', 1e-16,'RelTol', 1e-16,'NonNegative',1);
[T1,Y1] = ode15s(@ODEs2solve_MM1,[0 TimeLen],y0,options,piny,NSpecies);

%% Plot circulating drug
% Convert mol/mL to mol/L
Y1 = Y1 .* 1e3;
YT = Y1(:,1)+Y1(:,3);

figure; hold all;
% Blood unbound
plot(T1,log10(Y1(:,1)),'k');
% Blood bound
plot(T1,log10(Y1(:,3)),'r');
% Total bound + unbound
plot(T1,log10(YT(:,1)),'m');

legend('blood unbound','blood bound','blood total');
xlabel('Time(day)');
xlim([0 TimeLen]);
set(gca,'xtick',0:60*24:60*80);
set(gca,'xticklabel',0:1:60*80/24);
title('Drug conc over time mol/L');

% Make an array that summarizes the simulation ([time in {minutes-4},species]) 
Tvec = [5:TimeLen];
tmat = zeros(length(Tvec),NSpecies);
for nn = 1 : NSpecies
    tmat(:,nn)= interp1(T1,Y1(:,nn),Tvec);
end

% array ([time,summed drug (unbound and various bound) by distance]
for nn = 1 : length(Tvec)
    tsum(nn,:)=tmat(nn,NP+1:NP+nDB)...
        +FF*tmat(nn,NP+nDB*2+1:NP+nDB*3)...
        +tmat(nn,NP+nDB*3+1:NP+nDB*4)...
        +tmat(nn,NP+nDB*4+1:NP+nDB*5)...
        +tmat(nn,NP+nDB*6+1:NP+nDB*7);
end

cost = 0;
  