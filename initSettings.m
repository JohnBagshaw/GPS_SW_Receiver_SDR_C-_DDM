function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".  
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess = 36; %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels  = 8;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
settings.skipNumberOfBytes = 0;

% Code frequency basis
settings.codeFreqBasis     = 1.023e6;      %[Hz]

% Define number of chips in a code period
settings.codeLength         = 1023;

%% Raw signal file specific settings ===============================
cfID                       = fopen('datasetList.txt', 'r');
fileName                   = fscanf(cfID, '%c');
filePath                   = '';
settings.fileName          = [filePath, fileName];
[settings.fileID, message] = fopen(settings.fileName, 'rb');

if settings.fileID ~= -1
    fprintf('Probing data (%s)...\n', settings.fileName);
else
    error(['Error reading file: %s \nCheck contents of datasetList.txt file\n Check if'...
    ' specified file is present in relative path %s \n Check if file name specified has'...
    ' correct spelling\n%s'], fullFilePath, filePath, message);
end

% Specify acquisitin milliseconds
%settings.acqMS           = 1000;
settings.acqMS           = 2;
%if contains(settings.fileName, 'sample_nottochange_ch1_fileN')
if contains(settings.fileName, 'Woodbine_47a')
    
    settings.IF                   = 10e6;    %[Hz]
    settings.samplingFreq         = 40e6;    %[Hz]    
    settings.dataType             = 'int16'; % Data type used to store one sample
    settings.isBasebandSignal     = 1;
    
    settings.samplesPerCode       =  round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));
    settings.dataLengthBytes      = 3200000;
        
elseif contains(settings.fileName, 'GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin')
    
    settings.IF               = 4.1304e6;  %[Hz]
    settings.samplingFreq     = 16.3676e6; %[Hz]
    settings.dataType         = 'int8';    % Data type used to store one sample
    settings.isBasebandSignal = 0;
    
    settings.samplesPerCode    =  round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));
    settings.dataLengthBytes   = round(1145760000)/40; % MATLAB-safe length
    %settings.dataLengthBytes  =       1145760000;     % Actual length

elseif contains(settings.fileName, 'GPSdata-DiscreteComponents-fs38_192-if9_55.bin')
    
    settings.IF               = 9.548e6;   %[Hz]
    settings.samplingFreq     = 38.192e6;  %[Hz]
    settings.dataType         = 'int8';   % Data type used to store one sample
    settings.isBasebandSignal = 0;
    
    settings.samplesPerCode   =  round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));
    settings.dataLengthBytes  = round(1912602624)/40; % MATLAB-safe length
    %settings.dataLengthBytes =       1912602624;     % Actual length    
    
elseif contains(settings.fileName, 'GPSdata-DiscreteComponents-fs38_192-if9_55.bin')

    settings.IF               = (1590-1575.42)*1e6;   %[Hz]
    settings.samplingFreq     = 53e6;  %[Hz]
    settings.dataType         = 'ubit2';   % Data type used to store one sample
    settings.isBasebandSignal = 0;
    
    settings.samplesPerCode   =  round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));
    settings.dataLengthBytes  = settings.samplesPerCode*11; % MATLAB-safe length
    %settings.dataLengthBytes =       1912602624;     % Actual length    

else
    
    settings.IF               = (1590-1575.42)*1e6;   %[Hz]
    settings.samplingFreq     = 53e6;  %[Hz]
    settings.dataType         = 'ubit2';   % Data type used to store one sample
    settings.isBasebandSignal = 0;
    
    settings.samplesPerCode   =  round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));
    settings.dataLengthBytes  = settings.samplesPerCode*11; % MATLAB-safe length
    %settings.dataLengthBytes =       1912602624;     % Actual length    
    
    
end

if contains(settings.dataType, 'int8')
    settings.dataLength           = settings.dataLengthBytes / 1;
elseif contains(settings.dataType, 'int16')
    settings.dataLength           = settings.dataLengthBytes / 2;
elseif contains(settings.dataType, 'int32')
    settings.dataLength           = settings.dataLengthBytes / 4;
elseif contains(settings.dataType, 'ubit2')
    settings.dataLength           = settings.dataLengthBytes / 1;
else
    settings.dataLength           = settings.dataLengthBytes*4;
end

%% Acquisition settings ===================================================
% Algorithm type for acquisition
settings.acquisitionType = 'normalAcquisition';
% settings.acquisitionType = 'weakAcquisition';
% settings.acquisitionType = 'halfBitAcquisition';

% Decide if modulation is required based on 
% what acquisition algortihm accepts and format of dataset
if isequal(settings.acquisitionType, 'normalAcquisition')
   
    settings.dataExtractLen       = (11 + settings.acqMS) * settings.samplesPerCode * 4;
    % normal and halfbit acquisition accept IF data
    if  settings.isBasebandSignal == 1
        settings.modulationRequired   = 1;
        settings.demodulationRequired = 0;
    else
        settings.modulationRequired   = 0;
        settings.demodulationRequired = 0;
    end
    
elseif isequal(settings.acquisitionType, 'halfBitAcquisition')

    settings.dataExtractLen       =  settings.acqMS * settings.samplesPerCode;
    % normal and halfbit acquisition accept IF data
    if  settings.isBasebandSignal == 1
        settings.modulationRequired   = 1;
        settings.demodulationRequired = 0;
    else
        settings.modulationRequired   = 0;
        settings.demodulationRequired = 0;
    end
    
elseif isequal(settings.acquisitionType, 'weakAcquisition')
    
    settings.dataExtractLen       = settings.dataLength;
    % weak acquisition accept baseband data
    if  settings.isBasebandSignal == 1
        settings.modulationRequired   = 0;
        settings.demodulationRequired = 0;
    else
        settings.modulationRequired   = 0;
        settings.demodulationRequired = 1;
    end
    
end
% 
% if settings.dataExtractLen > settings.dataLength
%     warning('acqMS exceeds the duration of dataset.\nLimiting to maximum duration of dataset: %d',...
%         settings.dataLength);
%     prevDataExtractLen = settings.dataExtractLen;
%     prevAcqMS = settings.acqMS;
%     settings.dataExtractLen = settings.dataLength;
%     len = settings.dataExtractLen / settings.samplesPerCode / 2;
%     if prevAcqMS > floor(settings.dataLength/settings.samplesPerCode/2)
%         settings.acqMS = floor(len);
%     end
%     warning('acqMS reduced from %d to %d, dataExtractLen reduced from %d to %d', ...
%         prevAcqMS, settings.acqMS, prevDataExtractLen, settings.dataExtractLen);
% end


% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
settings.acqSatelliteList   = 1:32;         %[PRN numbers]
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 50;           %[kHz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2.5; % [MUST] Don't change it.

%% Tracking loops settings ================================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;       %[Hz]
settings.dllCorrelatorSpacing    = 0.5;     %[chips]

% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 25;      %[Hz]

%% Navigation solution settings ===========================================

% Period for calculating pseudoranges and position
settings.navSolPeriod       = 500;          %[ms]

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 10;           %[degrees 0 - 90]
% Enable/dissable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off
                                                          % 1 - On

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off
                                                          % 1 - On

%% Constants ==============================================================

settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time
