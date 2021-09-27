%--------------------------------------------------------------------------
%                           SDRJohn v1.0
% 
% Copyright (C) John Bagshaw
% Written by John Bagshaw
%--------------------------------------------------------------------------
%Not to be distributed to non-members of GNSS-lab, ESSE Department, LSE
%York University, Canada.
%--------------------------------------------------------------------------
%Script plots DDMs for acquisition results from C++ version of the software
%receiver.
%----------------------------------------------------------------------
%% Clean up the environment first =========================================
 clear; close all;
% For C++ 
fileName = ['NTLab_Bands_GPS_GLONASS_L12'];
load(['../GNSS_SDR_C++/',fileName,'_ddm_mat.mat']);

% This script is for plotting of DDM processing results.
% It takes acqResults by running DDMProcessing and plots
% the obtained correlation map and delay doppler map for
% all acquired PRNS.
if acquisitionType == 0
    acquisitionType = 'normalAcquisition';
else
    acquisitionType = 'weakAcquisition';
end

if isequal(acquisitionType, 'weakAcquisition')
    Nfd        = acqSearchBand / 0.5 + 1;
    Sblock     = floor((samplesPerCode * PIT) / Nfd);
    frqBins = IF - (acqSearchBand/2) * 1000 + 0.5e3 .* (0:(Nfd-1));
    
elseif isequal(acquisitionType, 'normalAcquisition')
    numberOfFrqBins = round(acqSearchBand * 2) + 1;
    
    for frqBinIndex = 1:numberOfFrqBins
        frqBins(frqBinIndex) = IF - ...
            (acqSearchBand/2) * 1000 + ...
            0.5e3 * (frqBinIndex - 1);
    end
end
ts = 1 / samplingFreq;

DelayDopplerMap = reshape(DelayDopplerMap, length(DelayDopplerMap)/32, 32);
DelayDopplerMap = {};

for i=1:32
    DelayDopplerMap{end+1,1} = reshape( ...
        samplesPerCode, numberOfFrqBins)';
end


% Get plot numbers
plotOffset = 1023;
CorrPlotIndex = 1 + plotOffset;
DDMPlotIndex = 2 + plotOffset;

% Option to change the set of PRNs
% These must be the acquired PRNs only.
useAllAcqPrns = 1;

acqSatelliteList = 1:32;
if useAllAcqPrns == 1
    acquiredPRNSubSet = acqSatelliteList;
else
    acquiredPRNSubSet = [4, 8];
end

frequencyBinIdx = frequencyBinIndex;
codePhas = codePhase;
for PRN = acquiredPRNSubSet
    
    % Select PRNs for which the threshold 
    % was exceeded for peak metric i.e. select
    % all acquired PRNs.
    
    if carrFreq(PRN) > 0

            results = DelayDopplerMap{PRN};
            
            
            figure(CorrPlotIndex);
            subplot(2,1,1)
            mesh(results)
            xlabel('Code phase [chips]')
            ylabel('Frequency [MHz]')
            title(['correlation Map for PRN ', num2str(PRN)]);
            
            subplot(2,1,2)
            imagesc(results)
            colorbar
            
            figure(DDMPlotIndex);
            subplot(2,1,1)
            
            dopplerAxis = frequencyBinIndex-20:frequencyBinIndex+20;
            
            if isequal(acquisitionType, 'weakAcquisition')
                dopplerAxis(find(dopplerAxis>Nfd)) = [];
            else
                dopplerAxis(find(dopplerAxis>numberOfFrqBins)) = [];
            end
            
            delayAxis = codePhase-50:codePhase+100;
            if isequal(acquisitionType, 'weakAcquisition')
                delayAxis(find(delayAxis > Sblock)) = [];
            else
                delayAxis(find(delayAxis > samplesPerCode)) = [];
            end
            
            DopplerY = frqBins(dopplerAxis)-carrFreq(PRN);
            TotalDelaysec = ((1:samplesPerCode)-codePhase)*ts*1e9;
            DelayX = TotalDelaysec(delayAxis);
            [Xax,Yax]=meshgrid(DelayX,DopplerY);
            DDM = results(dopplerAxis,delayAxis);
            
            
            mesh(Xax,Yax,DDM)
            xlabel('Delay axis [ns]')
            ylabel('Doppler axis [Hz]')
            title(['Delay Doppler Map for PRN ', num2str(PRN)]);
            box on
            
            subplot(2,1,2)
            imagesc(DelayX,DopplerY,DDM)
            colorbar
            xlabel('Delay axis [ns]')
            ylabel('Doppler axis [Hz]')
            title(['Delay Doppler Map for PRN ', num2str(PRN)]);
            
            fprintf('Plotted DDM for PRN %d\n', PRN);
            
    end 
end