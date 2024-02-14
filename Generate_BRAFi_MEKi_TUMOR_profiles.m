%% Master script for generating simulated BRAFi/MEKi tumor drug penetration model
%% Miles Miller, Thomas Ng

%% Common TUMOR PARAMETERS
input.RCAP    = 15e-4;                                                        %  capillary radius 15um
input.q       = log(2)./15;                                                   % absorption rate
input.deltaX  = 10 .* 10^-4;                                                  % integrating distance within the tumor in cm 

input.nDB = (500e-4) ./ input.deltaX;                                               % Rmax max radius of model, 500um

input.ATP_conc = 1e-3/1e3;                                                    % Intracellular ATP concentration

input.k_cat_ATP  = 0.66 * 60;                                                 % catabolism of ATP-bound moieties
input.k_unbind_ATP = 6.2 * 60;                                                % ATP-BRAF dissociation

input.TimeLen = 24*60*15;                                                     % Length of simulation (min)
%% Dabrafenib dependent parameters
input.DOSEmgkg  = 150;                                                        % Dose administered in mg
input.BIO       = 1;                                                          % Oral bioavailability (fraction)
input.Mass      = 1;                                                          % Mass of subject (Set as 1 for human, mouse mass in kg)
input.MW        = 520;                                                        % MW in g/mol

input.V1        = 70000*1.5;                                                  % volume of distribution (mL)  

input.k_bind_ATP = 0.11 * 1e6 * 1e3 * 60;                                     % ATP-BRAF association


input.DoseInt = 12*60;                                                        % Dosing interval (min)

input.k_bind = (1)*1.1e7 * 1e3;                                               % drug target association
input.k_unbind = (1)*0.006;                                                   % drug target dissociation

input.k_el1   = 0.012;                                                        % fast phase plasma elimination (1/min)
input.k_el2   = log(2) ./ (60*8);                                             % slow phase plasma elimination (1/min))
input.k_elA   = 0.3;                                                          % slow phase contribution (fraction)

input.kPerm1  = 3e-4 * 60 ;                                                   % permeability of tumor vasculature to unbound drug

input.Deff      = 64e-7 .* 60;                                                % effective diffusion of unbound drug

input.k_uptake  = 10*60;                                                      % partitioning into neutral lipid min-1
input.k_release = input.k_uptake ./ (10.^(3.65));                                   % partitioning into aqueous (lipid release)

input.DeffA     = 1e-7 .* 60;                                                 % Effective diffusion of albumin bound drug
input.kPerm1A   = 1e-7 * 60 ;                                                 % permeability of tumor vasculature to albumin
input.konALB    = (10^6 * 1E3 * 60) * (5e-4 * 1e-3);                          % albumin-drug association
input.koffALB   = input.konALB / (5e-4 * 1e-3) * (2e-6 * 1E-3 * 1 *1);              % albumin drug dissociation, calculated from KD and known drug plasma protein binding

input.k_uptakeALB = 0.001;                                                    % cellular drug-albumin uptake
input.FF = 0.07 + 0.3 * 0.01;                                                 % fraction of neutral lipid
input.k_catabALB = log(2) ./ (180);                                           % drug release from cellular albumin-drug complex

%%
[DAB_tmat, DAB_tsum, DAB_tfrac, ~, ~, ~] = simulate_TUMOR(input);

%% Trametinib dependent parameters
input.DOSEmgkg  = 2;                                                        % Dose administered in mg
input.BIO       = 0.7;                                                          % Oral bioavailability (fraction)
input.Mass      = 1;                                                          % Mass of subject (Set as 1 for human, mouse mass in kg)
input.MW        = 615;                                                        % MW in g/mol

input.V1        = 213087*0.5;                                                   % volume of distribution (mL)  

input.k_bind_ATP = 0;                                                        % ATP-BRAF association
                                                
input.DoseInt = 24*60;                                                        % Dosing interval (min)

input.k_bind = (1)*3.43E+05 * 1e3;                                               % drug target association
input.k_unbind = (1)*1.2E-4;                                                   % drug target dissociation

input.k_el1   = log(2) ./ (60*33);                                                        % fast phase plasma elimination (1/min)
input.k_el2   = log(2) ./ (60*4);                                             % slow phase plasma elimination (1/min))
input.k_elA   = 0.5;                                                          % slow phase contribution (fraction)

input.kPerm1  = 3e-4 * 60 ;                                                   % permeability of tumor vasculature to unbound drug

input.Deff      = 64e-7 .* 60;                                                % effective diffusion of unbound drug

input.k_uptake  = 10*60;                                                      % partitioning into neutral lipid min-1
input.k_release = input.k_uptake ./ (10.^(3.4));                                   % partitioning into aqueous (lipid release)

input.DeffA     = 1e-7 .* 60;                                                 % Effective diffusion of albumin bound drug
input.kPerm1A   = 1e-7 * 60 ;                                                 % permeability of tumor vasculature to albumin
input.konALB    = (10^6 * 1E3 * 60) * (5e-4 * 1e-3);                          % albumin-drug association
input.koffALB   = input.konALB / (5e-4 * 1e-3) * (2e-6 * 1E-3 * 6.25 *1);              % albumin drug dissociation, calculated from KD and known drug plasma protein binding

input.k_uptakeALB = 0.001;                                                    % cellular drug-albumin uptake
input.FF = 0.07 + 0.3 * 0.01;                                                 % fraction of neutral lipid
input.k_catabALB = log(2) ./ (180);                                           % drug release from cellular albumin-drug complex

%%
[TRAM_tmat, TRAM_tsum, TRAM_tfrac, ~, ~, ~] = simulate_TUMOR(input);

%% Encorafenib dependent parameters
input.DOSEmgkg  = 450;                                                        % Dose administered in mg
input.BIO       = 0.85;                                                          % Oral bioavailability (fraction)
input.Mass      = 1;                                                          % Mass of subject (Set as 1 for human, mouse mass in kg)
input.MW        = 540;                                                        % MW in g/mol

input.V1        = 163300*0.4;                                                   % volume of distribution (mL)  

input.k_bind_ATP = 0.11 * 1e6 * 1e3 * 60;                                                        % ATP-BRAF association
                                                
input.DoseInt = 24*60;                                                        % Dosing interval (min)

input.k_bind = 1.1E7;                                               % drug target association
input.k_unbind = 0.00039;                                                   % drug target dissociation

input.k_el1   = log(2) ./ (60);                                                        % fast phase plasma elimination (1/min)
input.k_el2   = log(2) ./ (60*3.5);                                             % slow phase plasma elimination (1/min))
input.k_elA   = 0.4;                                                          % slow phase contribution (fraction)

input.kPerm1  = 3e-4 * 60 ;                                                   % permeability of tumor vasculature to unbound drug

input.Deff      = 64e-7 .* 60;                                                % effective diffusion of unbound drug

input.k_uptake  = 10*60;                                                      % partitioning into neutral lipid min-1
input.k_release = input.k_uptake ./ (10.^(2.6));                                   % partitioning into aqueous (lipid release)

input.DeffA     = 1e-7 .* 60;                                                 % Effective diffusion of albumin bound drug
input.kPerm1A   = 1e-7 * 60 ;                                                 % permeability of tumor vasculature to albumin
input.konALB    = (10^6 * 1E3 * 60) * (5e-4 * 1e-3);                          % albumin-drug association
input.koffALB   = input.konALB / (5e-4 * 1e-3) * (2e-6 * 1E-3 * 40);              % albumin drug dissociation, calculated from KD and known drug plasma protein binding

input.k_uptakeALB = 0.001;                                                    % cellular drug-albumin uptake
input.FF = 0.07 + 0.3 * 0.01;                                                 % fraction of neutral lipid
input.k_catabALB = log(2) ./ (180);                                           % drug release from cellular albumin-drug complex

%%
[ENC_tmat, ENC_tsum, ENC_tfrac, ~, ~, ~] = simulate_TUMOR(input);

%% Binimetinib dependent parameters
input.DOSEmgkg  = 45;                                                        % Dose administered in mg
input.BIO       = 0.5;                                                          % Oral bioavailability (fraction)
input.Mass      = 1;                                                          % Mass of subject (Set as 1 for human, mouse mass in kg)
input.MW        = 441;                                                        % MW in g/mol

input.V1        = 91607*0.7;                                                   % volume of distribution (mL)  

input.k_bind_ATP = 0;                                                        % ATP-BRAF association
                                                
input.DoseInt = 12*60;                                                        % Dosing interval (min)

input.k_bind = 3.43E+05 * 1e3;                                               % drug target association
input.k_unbind = 1.57E-3;                                                   % drug target dissociation

input.k_el1   = 0.007;                                                        % fast phase plasma elimination (1/min)
input.k_el2   = log(2) ./ (60*3.5);                                             % slow phase plasma elimination (1/min))
input.k_elA   = 0.5;                                                          % slow phase contribution (fraction)

input.kPerm1  = 3e-4 * 60 ;                                                   % permeability of tumor vasculature to unbound drug

input.Deff      = 64e-7 .* 60;                                                % effective diffusion of unbound drug

input.k_uptake  = 10*60;                                                      % partitioning into neutral lipid min-1
input.k_release = input.k_uptake ./ (10.^(3.1));                                   % partitioning into aqueous (lipid release)

input.DeffA     = 1e-7 .* 60;                                                 % Effective diffusion of albumin bound drug
input.kPerm1A   = 1e-7 * 60 ;                                                 % permeability of tumor vasculature to albumin
input.konALB    = (10^6 * 1E3 * 60) * (5e-4 * 1e-3);                          % albumin-drug association
input.koffALB   = input.konALB / (5e-4 * 1e-3) * (2e-6 * 1E-3 * 10);              % albumin drug dissociation, calculated from KD and known drug plasma protein binding

input.k_uptakeALB = 0.001;                                                    % cellular drug-albumin uptake
input.FF = 0.07 + 0.3 * 0.01;                                                 % fraction of neutral lipid
input.k_catabALB = log(2) ./ (180);                                           % drug release from cellular albumin-drug complex

%%
[BINI_tmat, BINI_tsum, BINI_tfrac, ~, ~, ~] = simulate_TUMOR(input);

%% Vemurafenib dependent parameters
input.DOSEmgkg  = 960;                                                        % Dose administered in mg
input.BIO       = 0.6;                                                          % Oral bioavailability (fraction)
input.Mass      = 1;                                                          % Mass of subject (Set as 1 for human, mouse mass in kg)
input.MW        = 490;                                                        % MW in g/mol

input.V1        = 105548*2;                                                   % volume of distribution (mL)  

input.k_bind_ATP = 0.11 * 1e6 * 1e3 * 60;                                                        % ATP-BRAF association
                                                
input.DoseInt = 12*60;                                                        % Dosing interval (min)

input.k_bind = 1.1e7 *1E3;                                               % drug target association
input.k_unbind =0.03;                                                   % drug target dissociation

input.k_el1   = 0.0002434;                                                        % fast phase plasma elimination (1/min)
input.k_el2   = log(2) ./ (60*57);                                             % slow phase plasma elimination (1/min))
input.k_elA   = 0.5;                                                          % slow phase contribution (fraction)

input.kPerm1  = 3e-4 * 60 ;                                                   % permeability of tumor vasculature to unbound drug

input.Deff      = 64e-7 .* 60;                                                % effective diffusion of unbound drug

input.k_uptake  = 10*60;                                                      % partitioning into neutral lipid min-1
input.k_release = input.k_uptake ./ (10.^(5));                                   % partitioning into aqueous (lipid release)

input.DeffA     = 1e-7 .* 60;                                                 % Effective diffusion of albumin bound drug
input.kPerm1A   = 1e-7 * 60 ;                                                 % permeability of tumor vasculature to albumin
input.konALB    = (10^6 * 1E3 * 60) * (5e-4 * 1e-3);                          % albumin-drug association
input.koffALB   = input.konALB / (5e-4 * 1e-3) * (2e-6 * 1E-3 * 0.3);              % albumin drug dissociation, calculated from KD and known drug plasma protein binding

input.k_uptakeALB = 0.001;                                                    % cellular drug-albumin uptake
input.FF = 0.07 + 0.3 * 0.01;                                                 % fraction of neutral lipid
input.k_catabALB = log(2) ./ (180);                                           % drug release from cellular albumin-drug complex

%%
[VEM_tmat, VEM_tsum, VEM_tfrac, ~, ~, ~] = simulate_TUMOR(input);

%% Cobimetinib dependent parameters
input.DOSEmgkg  = 60;                                                        % Dose administered in mg
input.BIO       = 0.5;                                                          % Oral bioavailability (fraction)
input.Mass      = 1;                                                          % Mass of subject (Set as 1 for human, mouse mass in kg)
input.MW        = 531;                                                        % MW in g/mol

input.V1        = 802560*0.1;                                                   % volume of distribution (mL)  

input.k_bind_ATP = 0;                                                        % ATP-BRAF association
                                                
input.DoseInt = 24*60;                                                        % Dosing interval (min)

input.k_bind = 3.43E+05 * 1e3;                                               % drug target association
input.k_unbind =1.54E-3;                                                   % drug target dissociation

input.k_el1   = 0.002;                                                        % fast phase plasma elimination (1/min)
input.k_el2   = log(2) ./ (60*44);                                             % slow phase plasma elimination (1/min))
input.k_elA   = 0.9;                                                          % slow phase contribution (fraction)

input.kPerm1  = 3e-4 * 60 ;                                                   % permeability of tumor vasculature to unbound drug

input.Deff      = 64e-7 .* 60;                                                % effective diffusion of unbound drug

input.k_uptake  = 10*60;                                                      % partitioning into neutral lipid min-1
input.k_release = input.k_uptake ./ (10.^(3.9));                                   % partitioning into aqueous (lipid release)

input.DeffA     = 1e-7 .* 60;                                                 % Effective diffusion of albumin bound drug
input.kPerm1A   = 1e-7 * 60 ;                                                 % permeability of tumor vasculature to albumin
input.konALB    = (10^6 * 1E3 * 60) * (5e-4 * 1e-3);                          % albumin-drug association
input.koffALB   = input.konALB / (5e-4 * 1e-3) * (2e-6 * 1E-3 * 16.67);              % albumin drug dissociation, calculated from KD and known drug plasma protein binding

input.k_uptakeALB = 0.001;                                                    % cellular drug-albumin uptake
input.FF = 0.07 + 0.3 * 0.01;                                                 % fraction of neutral lipid
input.k_catabALB = log(2) ./ (180);                                           % drug release from cellular albumin-drug complex

%%
[COBI_tmat, COBI_tsum, COBI_tfrac, ~, ~, ~] = simulate_TUMOR(input);

%% Save data

save(fullfile('OUTPUT_DRUG_CONC.mat'),'TRAM_tsum', 'DAB_tsum', 'ENC_tsum', 'BINI_tsum','COBI_tsum', 'VEM_tsum' ,'TRAM_tmat', 'DAB_tmat', 'ENC_tmat', 'BINI_tmat','COBI_tmat', 'VEM_tmat', ...
    'ENC_tfrac', 'TRAM_tfrac', 'COBI_tfrac', 'BINI_tfrac', 'VEM_tfrac', 'DAB_tfrac');