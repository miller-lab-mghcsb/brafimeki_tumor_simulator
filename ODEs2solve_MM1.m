function dydt = ODEs2solve_MM1(t,y,pin,NSpecies2);

global SHIFT;
global DOSEDELAY DOSESTD nDB NP NCcount;
global BSITES BSITES_T;
global NSpecies;
global eps_d eps_a;

%%
%{
Input:
    t: list of timepoints
    y: list of concentration at each timepoint
    pin: parameter values
    NSpecies2: number of compartments/species to solve
Output:
    dydt: rate at which drug (CONCENTRATION) changes in each compartment
%}

% Initialize parameter values (see other doc for explanation)
for initializeparam = 1

q       = pin(1); 
V1      = pin(2);
k_el1   = pin(3);

kPerm1  = pin(4);
SAtumor = pin(5);

Vt  = pin(6);
Vs  = pin(7);
RCAP = pin(8);

DOSEmgkg  = pin(9);
MOUSEmass = pin(10);
DOSE      = pin(11);

Deff      = pin(12);
deltaX    = pin(13);

k_uptake  = pin(14);
k_release = pin(15);

KpT       = pin(16);

k_el2     = pin(17);
Kp_norm   = pin(18);
k_elA     = pin(19);
PermAnorm = pin(20);
Kp_ALBnorm = pin(21);
DeffA     = pin(22);
kPerm1A   = pin(23);
konALB    = pin(24);
koffALB   = pin(25);
k_uptakeALB = pin(26);
FF        = pin(27);

k_catabALB = pin(28);
k_bind     = pin(29);
k_unbind   = pin(30);

k_bind_ATP = pin(31);
k_unbind_ATP=pin(32);
ATP_conc   = pin(33);
k_cat_ATP  = pin(34);
TimeLen    = pin(35);
DoseInt    = pin(36);

end

dydt = zeros(NSpecies,1);

%% Compartment 1: drug in central compartment (plasma)

% Stomach absorption
dydt(1) = (q .* y(2) + k_el2 .* y(4)) ./ V1;

%one component elimination
dydt(1) = dydt(1) - k_el1 .* y(1);

% BOUND central drug
dydt(1) = dydt(1) - konALB .* y(1) + koffALB .* y(3);  

% BOUND central drug
dydt(3) = konALB .* y(1) - koffALB .* y(3);  %%Could make this a pseudo-equilibrium, unless off rate is slow (tight binding)

% BOUND clearance % same as unbound %% assumes rapid equilibrium 
dydt(3) = dydt(3) - k_el1 .* y(3);

% stomach
dydt(2) = - q .* y(2);
% If modeling repeat oral dosing...
for doseIndy = [60:DoseInt:TimeLen+60]
    dydt(2) = dydt(2) + DOSE .* normpdf(t,doseIndy,30);
end

% stomach slow
dydt(4) = - k_el2 .* y(4);
% If modeling repeat oral dosing...
for doseIndy = [60:DoseInt:TimeLen+60]
    dydt(4) = dydt(4) + k_elA * DOSE .* normpdf(t,doseIndy,30);
end


%% Compartment 2: perivascular tumor interstitium
NP = 4;
index = 1;

%method of lines
%reference: http://www.scholarpedia.org/article/Method_of_lines

% boundary condition no flux far from vessel
dydt(index+NP) = (Deff/deltaX.^2) * 6 * (y(index+NP+1) - y(index+NP));
dydt(index+NP+nDB*3) = (DeffA/deltaX.^2) * 6 * (y(index+NP+nDB*3+1) - y(index+NP+nDB*3));

% % Interstitial transport (method of lines)
for index = 2 : nDB - 1
    % free drug
    dydt(index+NP) = (Deff/deltaX.^2) * (y(index+1+NP) -2*y(index+NP) + y(index-1+NP));
    % dydt(index+NP) = dydt(index+NP) + (Deff ./ (2 * deltaX.^2 * (index-1)) * (y(index+1+NP)-y(index-1+NP)));
    dydt(index+NP) = dydt(index+NP) + (Deff ./ (deltaX.^2 * (index-1)) * (y(index+1+NP)-y(index-1+NP)));
          
    %ALB bound drug
    dydt(index+NP+nDB*3) = (DeffA/deltaX.^2) * (y(index+1+NP+nDB*3) -2*y(index+NP+nDB*3) + y(index-1+NP+nDB*3));
    %dydt(index+NP+nDB*3) = dydt(index+NP+nDB*3) + (DeffA ./ (2 * deltaX.^2 * (index-1)) * (y(index+1+NP+nDB*3)-y(index-1+NP+nDB*3)));
    dydt(index+NP+nDB*3) = dydt(index+NP+nDB*3) + (DeffA ./ (deltaX.^2 * (index-1)) * (y(index+1+NP+nDB*3)-y(index-1+NP+nDB*3)));

end

%transport [boundary condition] solving Neumann boundary dC/dr (at r=rcap) = P/D (Cp - Cc) 
dydt(nDB + NP) = Deff * ( y(nDB + NP - 1) - y(nDB + NP) + deltaX * kPerm1 * ( y(1) - y(nDB + NP)./eps_d ) ./ Deff ) * 2 / deltaX.^2;
% dydt(nDB + NP) = dydt(nDB + NP) + 2 * kPerm1 * ( y(1)-y(nDB + NP)./eps_d ) ./ Deff ./ (deltaX * (nDB-1));
dydt(nDB + NP) = dydt(nDB + NP) + 2 * kPerm1 * ( y(1)-y(nDB + NP)./eps_d ) ./ (deltaX * (nDB-1));

%BOUND transport into tumor insterstitium (SAME as unbound, same boundary,
%but different diffusion and permeability
dydt(nDB*4 + NP) = DeffA * ( y(nDB*4 + NP - 1) - y(nDB*4 + NP) + deltaX * kPerm1A * ( y(3) - y(nDB*4 + NP)./eps_a ) ./ DeffA ) * 2 / deltaX.^2;
% dydt(nDB*4 + NP) = dydt(nDB*4 + NP) + 2 * kPerm1A * ( y(3)-y(nDB*4 + NP)./eps_a ) ./ DeffA ./ (deltaX * (nDB-1));
dydt(nDB*4 + NP) = dydt(nDB*4 + NP) + 2 * kPerm1A * ( y(3)-y(nDB*4 + NP)./eps_a ) ./ (deltaX * (nDB-1));

for index = 1 : nDB
    
    % free drug [mol/mL]
    % dydt(index+NP) = dydt(index+NP) - FF * ( k_uptake * y(index+NP+nDB) * y(index+NP) ) + k_release * y(index+NP+nDB*2) * FF;
    dydt(index+NP) = dydt(index+NP) - FF * ( k_uptake * y(index+NP) ) + k_release * y(index+NP+nDB*2) * FF;
    dydt(index+NP) = dydt(index+NP) - konALB .* y(index+NP) ./ KpT + koffALB .* y(index+NP+nDB*3); 
    dydt(index+NP) = dydt(index+NP) + k_catabALB * y(index+NP+nDB*4) + koffALB * y(index+NP+nDB*4);
    dydt(index+NP) = dydt(index+NP) - k_bind * y(index+NP) * y(index+NP+nDB*5) + k_unbind * y(index+NP+nDB*6);
    
    % free lipid binding sites [mol/ mL] -- commented out b/c assuming
    % constant unsaturable lipid binding
    % dydt(index+NP+nDB) = - k_uptake * y(index+NP+nDB) * y(index+NP) + k_release * y(index+NP+nDB*2);
    
    % lipid bound drug [mol / mL]
    % dydt(index+NP+nDB*2) =  k_uptake * y(index+NP+nDB) * y(index+NP) - k_release * y(index+NP+nDB*2);
    dydt(index+NP+nDB*2) =  k_uptake * y(index+NP) - k_release * y(index+NP+nDB*2);
    
    % ALB bound drug
    dydt(index+NP+nDB*3) = dydt(index+NP+nDB*3) + konALB .* y(index+NP) ./ KpT - koffALB .* y(index+NP+nDB*3);
    dydt(index+NP+nDB*3) = dydt(index+NP+nDB*3) - k_uptakeALB * y(index+NP+nDB*3);

    % cellular ALB-bound drug
    dydt(index+NP+nDB*4) = k_uptakeALB * y(index+NP+nDB*3) - k_catabALB * y(index+NP+nDB*4) - koffALB * y(index+NP+nDB*4);
% % %     


    % free target binding sites [mol/L]
    dydt(index+NP+nDB*5) = - k_bind * y(index+NP) * y(index+NP+nDB*5) + k_unbind * y(index+NP+nDB*6) ...
        - ATP_conc * k_bind_ATP  * y(index+NP+nDB*5) + k_unbind_ATP * y(index+NP+nDB*7) + k_cat_ATP *  dydt(index+NP+nDB*7);
    
    % drug-bound target [mol / L]
    dydt(index+NP+nDB*6) =  + k_bind * y(index+NP) * y(index+NP+nDB*5) - k_unbind * y(index+NP+nDB*6);
    
     % ATP-bound target [mol / L]
   
    dydt(index+NP+nDB*7) =   ATP_conc * k_bind_ATP  * y(index+NP+nDB*5) - k_unbind_ATP * y(index+NP+nDB*7) - k_cat_ATP *  dydt(index+NP+nDB*7);

    % Phosphoylated downstream
%      dydt(index+NP+nDB*8)= k_cat_ATP *  dydt(index+NP+nDB*7) - (k_cat_ATP/100) *  dydt(index+NP+nDB*8);
     %This second term k_cat_ATP/100 is a placeholder for
     %dephosphorylation...arbitrary

end
