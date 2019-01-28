clc;
clear all;
close all;

NUM_BS = 1;
NUM_BS_ANTS = 16;
NUM_UE = 4;
FFT_SIZE_UE = 1024;
FFT_SIZE_BS = NUM_UE*FFT_SIZE_UE;
CP_SIZE = FFT_SIZE_BS/16;
SAMP_RATE_BS = 80e6;
SAMP_RATE_UE = 20e6;
NUM_SYMS = 100;

% OFDM modulator and preamble generation
ofdmMod = comm.OFDMModulator('FFTLength',FFT_SIZE_BS,'CyclicPrefixLength',CP_SIZE,'NumGuardBandCarriers',[0;0],'InsertDCNull',true,'NumSymbols',NUM_SYMS);
ofdm_demod = comm.OFDMDemodulator('FFTLength',FFT_SIZE_BS,'CyclicPrefixLength',CP_SIZE,'NumGuardBandCarriers',[0;0],'RemoveDCCarrier',true,'NumSymbols',NUM_SYMS);

% Raised Cosine Transmit Filter
UPSAMPLING_FACTOR = SAMP_RATE_BS / RATE; 
DOWNSAMPLING_FACTOR = SAMP_RATE_UE / RATE; 
DECIMATION_FACTOR = SAMP_RATE_UE / RATE;
FILTER_SPAN = 6;           % Filter span in symbol durations
BETA = 0.1;         % Roll-off factor

rctFilt = comm.RaisedCosineTransmitFilter(...
    'Shape',                  'Normal', ...
    'RolloffFactor',          BETA, ...
    'FilterSpanInSymbols',    FILTER_SPAN, ...
    'OutputSamplesPerSymbol', UPSAMPLING_FACTOR);

% Raised Cosine Receive Filter

FILTER_SPAN = 6;           % Filter span in symbol durations
BETA = 0.1;         % Roll-off factor

rcrFilt = comm.RaisedCosineReceiveFilter(...
    'Shape',                  'Normal', ...
    'RolloffFactor',          BETA, ...
    'FilterSpanInSymbols',    FILTER_SPAN, ...
    'InputSamplesPerSymbol', DOWNSAMPLING_FACTOR, ...
    'DecimationFactor', DECIMATION_FACTOR );

