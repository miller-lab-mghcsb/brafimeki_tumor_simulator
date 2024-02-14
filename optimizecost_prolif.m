%% Combined model derivation
% Miles Miller, Thomas Ng; 2024
% This script integrates the drug concentration vs. time vs. tumor location
% matrices with cell proliferation output.

clear all
load('PROLIF_DANA.mat');        % Load in vitro cell proliferation curves by Schulz et al. 2021

%Load output file from "Generate_BRAFi_MEKi_TUMOR_profiles.m"
load('OUTPUT_DRUG_CONC.mat');   % Load drug concentration matrices

timeperiod = 1:size(DAB_tsum,1); % Define the timeperiod to be analyzed.
tumorloc   = 1:size(DAB_tsum,2); % Define location within the tumor to be analyzed

% Optimize BRAFi/MEKi


BRAF_TFRAC{1}.TFRAC = DAB_tsum(timeperiod, tumorloc);
BRAF_TFRAC{2}.TFRAC = ENC_tsum(timeperiod, tumorloc);
BRAF_TFRAC{3}.TFRAC = VEM_tsum(timeperiod, tumorloc);

MEK_TFRAC{1}.TFRAC = TRAM_tsum(timeperiod, tumorloc);
MEK_TFRAC{3}.TFRAC = COBI_tsum(timeperiod, tumorloc);
MEK_TFRAC{2}.TFRAC = BINI_tsum(timeperiod, tumorloc);

conc = [100 10 1 1/8 1/64 1/(8*64) 1/(8^4) 1/(8^5) 1/(8^100)]; %MEKi concentration
conc2= conc.*10; %BRAFi concentration

lconc = log(conc);
lconc2= log(conc2);

for i = 1:3
    for j = 1:3
        cur_BRAF = BRAF_TFRAC{i}.TFRAC;
        cur_MEK  = MEK_TFRAC{j}.TFRAC;
        cur_PROLIF=PROLIF{i,j}.V;

        %Extrapolate
        new_PROLIF = zeros(size(cur_PROLIF,1)+2, size(cur_PROLIF,2)+2);
        new_PROLIF(3:end,3:end) = cur_PROLIF;


        for k = 1:size(cur_PROLIF,1)
            new_PROLIF(2,k+2) = interp1(lconc2(3:4), cur_PROLIF(1:2,k), lconc2(2), 'linear', 'extrap');
        
            new_PROLIF(1,k+2) = interp1(lconc2(2:3), new_PROLIF(2:3,k+2), lconc2(1), 'linear', 'extrap');
        
        end

         for k = 1:size(new_PROLIF,2)
            new_PROLIF(k,2) = interp1(lconc(3:4), new_PROLIF(k,3:4), lconc(2), 'linear', 'extrap');
        
            new_PROLIF(k,1) = interp1(lconc(2:3), new_PROLIF(k,2:3), lconc(1), 'linear', 'extrap');
        
        end

        ind = find(new_PROLIF < 0);
        new_PROLIF(ind) = 0;

         ind = find(new_PROLIF > 100);
        new_PROLIF(ind) = 100;


        cur_BRAF = cur_BRAF.*1E6;
        cur_MEK  = cur_MEK.*1E6;

%         cur_BRAF(cur_BRAF >max(conc2)) = max(conc2);
%         cur_MEK(cur_MEK > max(conc))   = max(conc);

        cur_BRAF = cur_BRAF(:,:);
        cur_MEK  = cur_MEK(:,:);

        [X,Y] = meshgrid(conc, conc2');
        X = X'; %MEK
        Y = Y'; %BRAF

       % F = griddedInterpolant(X,Y, cur_PROLIF);

      
%             cur_BRAF(cur_BRAF == 0) = 1/(8^100);
%               cur_MEK(cur_MEK == 0) = 1/(8^100);
       

        Vq = interp2(Y, X, new_PROLIF, cur_BRAF, cur_MEK);

     % F = griddedInterpolant(conc, conc2', cur_PROLIF);
      
        Vq(Vq < 0) = 0;
        Vq(Vq > 100) = 100;
        metric{i,j}.Vq = Vq;

    end
end

meanMETRIC = zeros(3,3);
maxMETRIC  = zeros(3,3);
minMETRIC  = zeros(3,3);
STDMETRIC  = zeros(3,3);

for i = 1:3
    for j = 1:3
        curVq = metric{i,j}.Vq;
        meanMETRIC(i,j) = mean(curVq(:));
        maxMETRIC(i,j) = max(curVq(:));
        minMETRIC(i,j) = min(curVq(:));
        STDMETRIC(i,j) = std(curVq(:));
    end
end

% Row: BRAFi (Dab, Enc, Vem); Column:  MEKi (Tram, Cobi, Bini);
save('Integrated_Model_output.mat', 'metric', 'meanMETRIC', 'maxMETRIC', 'minMETRIC', 'STDMETRIC');

