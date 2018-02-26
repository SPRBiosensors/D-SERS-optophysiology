function [eventOpto, Scores, BarcodeStandard] = Optophysiology(Wave, Data,SERSDatab, varargin)    
% %--------------------------------------------------------------------------
% % Version 3.0 (Updated 2017 August 20)
% % Optophysiology(Wave, Data, SERSDatab, Options)
% %   Generate the OptoPhysiology curves, BarCodes of the Statdard and number
% %   of events for the probed metabolites.
% % 
% %   Inputs:     Wave                        Vector of Raman Shift.
% %               Data                        Matrix of spectra (Kinetic).
% %               SERSDatab                   Cell containing the name and
% %                                           the chemometric parameters for
% %                                           the sorting.
% % 
% %   Options:    'all'                       Analysis carried on all standards in
% %                                           SERSDatab.
% %               'Name'                      Chemometric only on 'Name'
% %                                           standard in SERSDatab.
% % 
% %   References:
% %           (1)  Eilers, P., H., C., A perfect smoother, Analytical Chemistry 75 (14), 3631 (2003).
% %           (2)  Zhang, Z. M., Chen, S., Liang, Y.-Z., Analyst, 2010, 135, 1138-1146
% %           (3)  Butler, H. J. et al, Nature Protocols, 2016, 11 (4), 664-687 
% %             
% % 
% %   Felix Lussier, Thibeault Brule, Université de Montréal, July 2017
%--------------------------------------------------------------------------
% OPEN
shift = Wave; serie = Data; echanti = 20; echant2 = 8; echantiRef = 20;
echant = echanti.*1024./(max(shift)-min(shift)); echant = round(echant);
nbspec = size(serie,2); sizwave = size(shift,2); 
%--------------------------------------------------------------------------
% STANDARDS BARCODES EXTRACTION
% Import Standard Spectrum and parameter from external .mat file
[Standard,shiftRef,nbref,RefName] = StandardImport(SERSDatab,varargin{:});

shiftRef = shiftRef(235:end,:); %Standard acquiered on a previous Raman microscope with a different resolution. Spectrum must be truncted.
DatabFits = zeros(nbref,size(shiftRef,1)); 
Serie = zeros(sizwave,nbspec); SerieFits = zeros(nbspec,sizwave);
mph = zeros(sizwave,nbspec); mphRef = zeros(size(shiftRef,1),nbref);
scores = zeros(nbspec,nbref); Scores = scores; eventOpto = cell(nbref,2);
Nbpics = cell2mat(SERSDatab(:,2)); Nbpicsref = cell2mat(SERSDatab(:,3));
BarcodeStandard = zeros(sizwave,nbref);

Lambda_PLS = 1e8; %Value of Lambda used in the penalized least square regression for the background (Ref 1-2)

datab = zeros(size(shiftRef,1),nbref);
Datab = datab;
diffDatab = zeros(size(shiftRef,1)-2, nbref);
diffSerie = zeros(sizwave - 2, nbspec);

% Preprocessing of the Standard spectra.
for h = 1:nbref;
    %Vector normalization can be added.
    Datab(:,h) = sgolayfilt(Standard(:,h),3,11); %Noise filter on Standard
    [~,DatabFits(h,:)] = airPLS(Datab(:,h)', Lambda_PLS, 2, 0.05); %Extraction of the background. Downloaded from Ref 2.
    mphRef(:,h) = 0.0055*median(Datab(:,h)) + DatabFits(h,:)'; %Curved Threshold line
    frtdiff = diff(Datab(:,h)); %First differentiation
    scnddiff = diff(sgolayfilt(frtdiff,3,11)); %Second diff.
    scnddiff(scnddiff > 0) = 0; %Extraction of ONLY the Raman peak, no shoulders
    diffDatab(:,h) = -scnddiff;
end    

% Preprocessinf of the Experiemental Spectra.
for i = 1:nbspec;
    %Vector Normlaization can be added.
    Serie(:,i) = serie(:,i);
    [~,SerieFits(i,:)] = airPLS(Serie(:,i)', Lambda_PLS, 2, 0.05); %Extraction of the background
    Serie(:,i) = sgolayfilt(Serie(:,i),3,11); %Noise filter
    mph(:,i) = 0.0055*median(Serie(:,i)) + SerieFits(i,:)'; %Curved Threshold line
    frtdiff = diff(Serie(:,i)); %First differentation
    scnddiff = diff(sgolayfilt(frtdiff,3,11)); %Second diff
    scnddiff(scnddiff > 0) = 0; %Extraction of ONLY the Raman peak, no shoulders
    diffSerie(:,i) = -scnddiff;
end

%End preprocessing
% -------------------------------------------------------------------------
%Loop for anayling all the experimental spectra:
    
for h = 1:nbref; nbpics = Nbpics(h,1); nbpicsref = Nbpicsref(h,1);
    for i = 1:nbspec;
    
    echantRef = echantiRef.*1024./(max(shiftRef(:,h)) - min(shiftRef(:,h))); 
    echantRef = round(echantRef);
    
    clear picsref lcsref pics lcs
         
    [picsref,lcsref]=findpeaks(diffDatab(:,h),'sortstr','descend',...
                                          'minpeakdistance',echantRef); %Localized Raman peaks in standards.
      
    picsref = picsref'; lcsref = lcsref';                                 
    if size(picsref,2)~= nbpicsref;
        Lcsref = horzcat(lcsref,zeros(1,nbpicsref)); Lcsref = Lcsref(:,1:nbpicsref);   
    else Lcsref = lcsref;
    end; locsref = Lcsref;
    

    for j = 1:size(locsref,2);
        if Datab(locsref(:,j),h) < mphRef(locsref(:,j),h); %Applied Threshold
            locsref(:,j) = 0;
        end
    end
    
    for z = 1:size(locsref,2);
        if locsref(1,z) ~= 0;
           idx = find(shiftRef(:,h) >= shiftRef(locsref(1,z),h) - echant2 &...
                 shiftRef(:,h) <= shiftRef(locsref(1,z),h) + echant2);
                 for k = 1:size(idx,1);
                     BarcodeStandard(idx(k),h) = max(Standard(:,h)); %Generation of the Standards Barcodes.
                 end
        end
    end 
    
        [pics,lcs] = findpeaks(diffSerie(:,i),'sortstr','descend',...
                                          'minpeakdistance',echant); %Localized Raman peaks in the experimental spectrum
                                      
        pics = pics'; lcs = lcs';
                                      
        if size(pics,2) ~= nbpics;
            Lcs = horzcat(lcs,zeros(1,nbpics)); Lcs = Lcs(:,1:nbpics);
        else Lcs = lcs;
        end; locs = Lcs; 
            
        for j = 1:size(locs,2);
            if Serie(locs(:,j),i) < mph(locs(:,j),i); %Applied Threshold
                locs(:,j) = 0;
            end
        end
%--------------------------------------------------------------------------
%Barcodes comparision (Standard VS Experiemental). Notes that it still
%within the Loop.

        ident = zeros(1,nbpics);
            for j = 1:nbpics;
                for l = 1:nbpicsref; %Convert both Barcodes from pixel to Raman Shift prior comparision. Allow Universal comparision.
                    if  locsref(1,l)~=0 && locs(1,j) ~= 0 &&...
                        shiftRef(locsref(1,l),h) - echant2 <= shift(locs(1,j)) &&...                            
                        shift(locs(1,j)) <= shiftRef(locsref(1,l),h) + echant2; % Check, one by one, if the bar is in the acceptability threshold (width of tha bar)
                        ident(1,j) = 1; % If the bar is in the acceptability range, a positive value (1) is obtained.
                    end
                end
            end
            scores(i,h) = sum(ident(1,:),2); % Final result correspond to the total number of positive bar per experimental barcode for every Standard.
    end
    
    for j = 1:size(scores,1);
        if scores(j,h) >= round(0.6*nbpicsref); % Apply a threshold on the number of positive bar required for a positive event, for every standard.
           Scores(j,h) = 1;
            else Scores(j,h) = 0;
        end
    end
end
    
for i = 1:nbref; %Display the name, and  the counts for the probed metabolites.
    eventOpto{i,1} = RefName{i,1}; 
    eventOpto{i,2} = sum(Scores(:,i)); 
end;
% CLOSE
%--------------------------------------------------------------------------




