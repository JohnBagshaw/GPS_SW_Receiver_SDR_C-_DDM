function acqResults = acquisition(longSignal, settings)
% load('longsignal_woodbine2_47_2ms.mat');
% load('settings_woodbine2_47_2ms.mat');
% longSignal = data;
% % %settings = settings;

%
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

signal0DC = longSignal - mean(longSignal); 

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable(settings);


%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList
    
%    PRN=31;

%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));
    
    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins

        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
        cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

        %--- "Remove carrier" from the signal -----------------------------
        acqResIfft1 = zeros(1, samplesPerCode);
        acqResIfft2 = zeros(1, samplesPerCode);
        
        if floor(settings.acqMS/2) == 0
            error('acqMS value (specified in initSettings.m)',...
                  'should be atleast 2 for normal acquisition to work.');
        end
        
        for ms = 1 : floor(settings.acqMS/2)
            
            msIdx1 = 2*(ms-1)+1;
            msIdx2 = 2*(ms-1)+2;
            
            % Create two 1msec vectors of data to correlate with and one with zero DC
            signal1 = longSignal(1 : samplesPerCode);
            signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);
            
            I1      = sinCarr .* signal1;
            Q1      = cosCarr .* signal1;
            I2      = sinCarr .* signal2;
            Q2      = cosCarr .* signal2;
            
            %--- Convert the baseband signal to frequency domain --------------
            IQfreqDom1 = fft(I1 + 1i*Q1);
            IQfreqDom2 = fft(I2 + 1i*Q2);
            
            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
            convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
            
            %--- Perform inverse DFT and store correlation results ------------
            acqResIfft1 = acqResIfft1 + ifft(convCodeIQ1);
            acqResIfft2 = acqResIfft2 + ifft(convCodeIQ2);
        end
        
        acqRes1 = abs(acqResIfft1) .^ 2;
        acqRes2 = abs(acqResIfft2) .^ 2;

        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
        else
            results(frqBinIndex, :) = acqRes2;
        end

    end % frqBinIndex
    
    % save results only when all frequency
    % search bins are processed
    acqResults.DelayDopplerMap(PRN) = {results};

    %% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [peakSize frequencyBinIndex] = max(max(results, [], 2));
    acqResults.frequencyBinIndex(PRN)=frequencyBinIndex;
    
    %--- Find code phase of the same correlation peak ---------------------
    [peakSize codePhase] = max(max(results));
    acqResults.codePhase(PRN)=codePhase;
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
        
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Generate 10msec long C/A codes sequence for given PRN --------
        caCode = generateCAcode(PRN);
        
        codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                           
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
    
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier = ...
            signal0DC(codePhase:(codePhase + 10*samplesPerCode-1)) ...
            .* longCaCode;
        
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        fftxc = abs(fft(xCarrier, fftNumPts)); 
        
        uniqFftPts = ceil((fftNumPts + 1) / 2);
        [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
        
        fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
        
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex);
        acqResults.codePhase(PRN) = codePhase;
    
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
%     figure(1023)
%     subplot(2,1,1)
%     mesh(results)
%     ylabel('Frequency')
%     xlabel('Code phase [chips]')
%     title('Delay Doppler Map')
%     subplot(2,1,2)
%     imagesc(results)
%     colorbar
    
end    % for PRN = satelliteList
%=== Acquisition is over ==================================================
% fprintf(')\n');
Results_DDM=cell2mat(acqResults.DelayDopplerMap(20));
[a,b]=size(Results_DDM);
Signal_MAX=max(max(Results_DDM));
Signal_Avg0=mean(mean(Results_DDM));
Signal_Avg1=(mean(mean(Results_DDM)).*(a*b)-Signal_MAX)./(a*b-1);
q=0;
qq=0;
for i=1:a
    for j=1:b
        if Results_DDM(i,j)<=(30/100)*Signal_MAX
            q=q+1;
        Results_DDM_Trailing_edg_noise(q,1)=Results_DDM(i,j);
        elseif Results_DDM(i,j)>(30/100)*Signal_MAX
            qq=qq+1;
        Results_DDM_Trailing_edg_signal(qq,1)=Results_DDM(i,j);    
        end
    end
    
end
Signal_Avg2=mean(Results_DDM_Trailing_edg_noise);
Signal_30=mean(Results_DDM_Trailing_edg_signal);
SNR_1=Signal_MAX./Signal_Avg1;
SNR_2=Signal_MAX./Signal_Avg2;
SNR_3=Signal_30./Signal_Avg2;
SNR_Baseband_wikipedia=(Signal_Avg0.^2)./((1./(a*b)).*(sum(sum((Results_DDM-Signal_Avg0).^2))));

% frequencyBinIndex=acqResults.frequencyBinIndex(12);
% codePhase=acqResults.codePhase(12);
% PRN=20;
%  figure(1023)
%     subplot(2,1,1)
%     mesh(Results_DDM(:,:))
%     xlabel('Code phase [chips]')
%     ylabel('Frequency [MHz]')
%     title('correlation Map')
%     
%     subplot(2,1,2)
%     imagesc(Results_DDM(:,:))
%     colorbar
%     
%     figddm = figure(1024);
%     subplot(2,1,1)
% 
%     dopplerAxis = frequencyBinIndex-10:frequencyBinIndex+10;   %-20 +20
%     delayAxis = codePhase-50:codePhase+100;                    %-50 +100
% 
%     DopplerY = frqBins(dopplerAxis)-acqResults.carrFreq(PRN);
%     TotalDelaysec = ((1:samplesPerCode)-codePhase)*ts*1e9;
%     DelayX = TotalDelaysec(delayAxis);
%     [Xax,Yax]=meshgrid(DelayX,DopplerY);
%     DDM = Results_DDM(dopplerAxis,delayAxis);
% 
%     
%     mesh(Xax,Yax,DDM) 
%     xlabel('Delay axis [ns]')
%     ylabel('Doppler axis [Hz]')
%     title('Delay Doppler Map')
%     box on
%     
%     subplot(2,1,2)
%     imagesc(DelayX,DopplerY,DDM)
%     colorbar
%     xlabel('Delay axis [ns]')
%     ylabel('Doppler axis [Hz]')
%     title('Delay Doppler Map')
% 
% %     filename = strsplit(settings.fileName,'/');
% %Savefilename = strsplit(filename{end},'.bin');
% %saveas(figddm,sprintf('%s.jpg',Savefilename{1}))
% %end    % for PRN = satelliteList
% %=== Acquisition is over ==================================================
% fprintf(')\n');