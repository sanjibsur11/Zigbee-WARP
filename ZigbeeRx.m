% Implementation of 802.15.4 PHY decoder
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments: This file contains the 802.15.4 PHY decoding function
% 
% Version: 0.0.2
% Last modified: 01/17/2014
% 
% Comments: Adding down sampling code for WARP board transmission. Modified
% the sensitivity of FM demodulator. Adding carrier sensing to remove
% signals not required. Removing clock recovery circuit parameters. Found a
% bug regarding Start Frame Delimeter identification. Identifier 0xA7 has
% to be identified by two consequtive 32 bit chips. In earlier parsing 0xA
% and 0x7 may be received discontiguously, which was not correct.
% 


function [BER, CER, syncFound] = ZigbeeRx(rxData_board_osamp)

global USESIM DEFAULT_ERROR_THRESHOLD CHIP_ENCODE_MAPPING ...
    CHIP_DECODE_MAPPING STATE_SYNC_SEARCH STATE_HAVE_SYNC ...
    STATE_HAVE_HEADER numBitsTotal numSymbolsTotal SIMSNR useFakeChannel...
    DEBUG_ON VERBOSE1 VERBOSE2 numZeroPreamble osamp padding_size;

%% Chip decoding map
            
% Chips decoding sequence, note that it is different from the original 
% encoding sequence, because of O-QPSK and half sine pulse filter at the
% transmitter side
CHIP_DECODE_MAPPING = ...
                [1618456172
                1309113062
                1826650030
                1724778362
                778887287
                2061946375
                2007919840
                125494990
                529027475
                838370585
                320833617
                422705285
                1368596360
                85537272
                139563807
                2021988657];

            

%% Downsampling
            
if USESIM
    rxData_board = rxData_board_osamp;
else
    % For WAPR Board transmission downsample Rx data by at least 2
    rxData_board = decimate(rxData_board_osamp, osamp);
end



%% Perform carrier sense to get the start of signal
idx = carrier_sense(rxData_board);
if idx ~= 0
    % Found valid signal
    % Now padding some zero at the beggining to help bootstrap the FM
    % demodulator and M&M logic
    rxData_board = [zeros(1, padding_size).'; rxData_board(idx:end)];
% end

if DEBUG_ON
    figure(98);
    plot(1:length(rxData_board), real(rxData_board), 'b-', ...
         1:length(rxData_board), imag(rxData_board), 'r-');
    xlabel('Sample count');
    ylabel('Rx samples');
    title('Received samples after carrier sensing');
end


%% FM demodulation

% Samples per symbol is 2 for FM demodulation and clock recovery
sps = 2;

% Convert complex samples to raw floating point samples using FM
% demodulation
sensitivity = (pi / 2) / sps;
% sensitivity = 1;
% Raw samples from FM demodulation contains DC offset
raw_samples_wDCO = fm_demod(rxData_board, 1/sensitivity);



%% Low pass filter

% Low pass the output of FM Demod to allow us to remove the DC offset
% resulting from the frequency offset
% Weight factor for single pole IIR filter to remove DC offset

if USESIM
    alpha = 0;
else
    alpha = 0;
end

% Removing the DC offset using frequency compensation
raw_samples = remove_dco(alpha, raw_samples_wDCO);
% raw_samples = raw_samples_wDCO;


%% Mueller and Müller (M&M) clock recovery logic
% Recovery related parameters
omega = sps;
if USESIM
	gain_mu = 0;
	mu = 0;
	omega_relative_limit = 0;
	% Critically damped
	gain_omega = 0;
  
else % If USESIM=0, Use WARP Board parameters
    gain_mu = 0;
    mu = 0;
    omega_relative_limit = 0;
    % Critically damped
    gain_omega = 0;
end


% Calling the Mueller and Müller (M&M) function to generate +/- float data
clock_recovered_samples = mueller_muller(omega, gain_omega, mu, gain_mu,...
    omega_relative_limit, raw_samples);



%% Binary slicer to convert +/- float values to binary 0/1
% Initialize decode related variables
skipsamples = 0;

sigin = clock_recovered_samples;
inbuf = (sigin > 0);



%% DSSS demodulation
% =========== Loop through all samples ============
ns = skipsamples;
zero_count = 0;
d_preamble_cnt = 0;
d_packet_byte = 0;
sync_found = 0;

fprintf('Starting the synchronization vector search process\n');  

while ns < length(inbuf) - 32
    
%% Start: Locate contiguous '0' samples and SFD
    
    
    if d_preamble_cnt == 0
        % Grab 32 chips
        d_shift_reg = inbuf(ns+1:ns+32);
        ns = ns + 1;
        % Decode chips
        [symbol error_count] = decode_chips(bi2de(d_shift_reg,...
                                                            'left-msb'));
        % Find whether matching with '0' or not
        if (symbol == 0) && (error_count <= DEFAULT_ERROR_THRESHOLD)
            zero_count = zero_count + 1;
            d_preamble_cnt = d_preamble_cnt + 1;
            fprintf('Found first 0 in the packet!\n');
            % ns is already rolled extra one time
            ns = ns - 1;
        else
            continue;
        end
    else
        % We found the first '0', now grab every 32 chips
        ns = ns + 32;
        
        % Check whether length is beyond the array
        if ns + 32 > length(inbuf)
            break;
        end
        
        d_shift_reg = inbuf(ns+1:ns+32);
        if d_packet_byte == 0
            [symbol error_count] = decode_chips(bi2de(d_shift_reg,...
                                                            'left-msb'));
            if error_count <= DEFAULT_ERROR_THRESHOLD
                if symbol == 0
                    zero_count = zero_count + 1;
                    d_preamble_cnt = d_preamble_cnt + 1;
                % Found the LSB of PACKET_SYNC_VECTOR
                elseif symbol == hex2dec('A')
                    fprintf('Found %d 0s in the packet!\n', d_preamble_cnt);
                    fprintf('Found first SFD nibble in packet!\n');
                    d_packet_byte = bitshift(symbol, 4);
                end
            else
                if VERBOSE2
                    fprintf('Error count is too high!\n');
                end
            end
        else
            [symbol error_count] = decode_chips(bi2de(d_shift_reg,...
                                                            'left-msb'));
            if error_count <= DEFAULT_ERROR_THRESHOLD
                if symbol == 7
                    d_packet_byte = d_packet_byte + symbol;
                    fprintf(['Found the complete synchronization header'...
                        ' 0x%x\n'], d_packet_byte);
                    sync_found = 1;
                else
                    if VERBOSE2
                        fprintf('Wrong byte in synchronization header\n');
                    end
                end
                break;    
            else
                if VERBOSE2
                    fprintf('Error count is too high!\n');
                end
            end
        end
    end
    
%% End: Locate SFD 
end



%% Synchronization header is found, now decode the rest of samples

if sync_found
    ns = ns + 32;
    
    if (ns + numBitsTotal*(32/4) > length(inbuf))
        inbuf = [inbuf zeros(1, ns + numBitsTotal*(32/4) - length(inbuf))];
    end
    
    % Every 4 bit is mapped to a 32 bit
    minbuf = inbuf(ns+1:ns+numBitsTotal*(32/4));
    inchips = reshape(minbuf, 32, numSymbolsTotal).';

    bit_error = 0;

    for ds = 1:numSymbolsTotal
        [symbol(ds) error_cnt(ds)] = decode_chips(bi2de(inchips(ds,:),...
                                                            'left-msb'));
    end
    
    symbol(symbol < 0) = 0;
    obits = de2bi(symbol, 'left-msb');
    for k = 1:numSymbolsTotal
        outbits(k,:) = [zeros(1, 4-length(obits(k,:))) obits(k,:)];
    end
    
    outbits = reshape(outbits.', 1, numBitsTotal).';
    inbits = load('databits.dat');
    mBER = sum(outbits ~= inbits)/length(inbits);
    mCER = sum(error_cnt)/(length(error_cnt)*32);
    
    fprintf('BER = %0.6g%%\n', mBER*100);
    fprintf('CER = %0.6g%%\n', mCER*100);
    BER = mBER;
    CER = mCER;
    syncFound = 1;
    
    if DEBUG_ON
        
        figure(101);
        stairs(inbits);
        axis([0 length(inbits) -0.5 1.5]);
        xlabel('Bit number');
        ylabel('Bit value');
        title('Input bits');
        
        figure(102);
        stairs(outbits);
        axis([0 length(outbits) -0.5 1.5]);
        xlabel('Bit number');
        ylabel('Bit value');
        title('Output bits');
        
        bit_error = xor(outbits, inbits);
        idx = (bit_error == 0);
        x = 1:length(bit_error);
        figure(103);
        plot(x(idx), bit_error(idx), 'g*', x(~idx), bit_error(~idx), 'ro');
        axis([0 length(bit_error) -0.5 1.5]);
        xlabel('Bit position');
        ylabel('Error');
        title('Error position');
        
    end
    
else
    fprintf('Synchronization header not found!\n');
    BER = 0;
    CER = 0;
    syncFound = 0;
end


fprintf('\n\n');

end