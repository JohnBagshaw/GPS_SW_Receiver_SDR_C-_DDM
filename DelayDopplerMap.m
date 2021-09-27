% This script is for plotting of DDM processing results.
% It takes acqResults by running DDMProcessing and plots
% the obtained correlation map and delay doppler map for
% all acquired PRNS.

if isequal(settings.acquisitionType, 'weakAcquisition')
    Nfd        = settings.acqSearchBand / 0.5 + 1;
    Sblock     = floor((samplesPerCode * settings.PIT) / Nfd);
    frqBins = settings.IF - (settings.acqSearchBand/2) * 1000 + 0.5e3 .* (0:(Nfd-1));
    
elseif isequal(settings.acquisitionType, 'normalAcquisition')
    numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;
    samplesPerCode = round(settings.samplingFreq / ...
        (settings.codeFreqBasis / settings.codeLength));
    for frqBinIndex = 1:numberOfFrqBins
        frqBins(frqBinIndex) = settings.IF - ...
            (settings.acqSearchBand/2) * 1000 + ...
            0.5e3 * (frqBinIndex - 1);
    end
end
ts = 1 / settings.samplingFreq;


% Get plot numbers
plotOffset = 1023;
CorrPlotIndex = 1 + plotOffset;
DDMPlotIndex = 2 + plotOffset;

% Option to change the set of PRNs
% These must be the acquired PRNs only.
useAllAcqPrns = 0;

if useAllAcqPrns == 1
    acquiredPRNSubSet = settings.acqSatelliteList;
else
    acquiredPRNSubSet = [4, 8];
end

for PRN = acquiredPRNSubSet
    
    % Select PRNs for which the threshold 
    % was exceeded for peak metric i.e. select
    % all acquired PRNs.
    
    if acqResults.carrFreq(PRN) > 0

            results = acqResults.DelayDopplerMap{PRN};
            frequencyBinIndex = acqResults.frequencyBinIndex(PRN);
            codePhase = acqResults.codePhase(PRN);
            
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
            dopplerAxis(find(dopplerAxis<1)) = [];
            if isequal(settings.acquisitionType, 'weakAcquisition')
                dopplerAxis(find(dopplerAxis>Nfd)) = [];
            else
                dopplerAxis(find(dopplerAxis>numberOfFrqBins)) = [];
            end
            
            delayAxis = codePhase-50:codePhase+100;
            if isequal(settings.acquisitionType, 'weakAcquisition')
                delayAxis(find(delayAxis > Sblock)) = [];
            else
                delayAxis(find(delayAxis > samplesPerCode)) = [];
            end
            
            DopplerY = frqBins(dopplerAxis)-acqResults.carrFreq(PRN);
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
            
            CorrPlotIndex = CorrPlotIndex + 2;
            DDMPlotIndex = DDMPlotIndex + 2;
            
    end 
end

