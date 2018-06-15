
clear

%% Generate BLE Channels

% Signal Parameters
Fs = 1e10;                                  % Sampling frequency
T = 1/Fs;                                   % Sampling period
L = 100000;                                 % Length of signal
t = (0:L-1)*T;                              % Time vector

% Transmit Power
Pt_dBm = 10;
Pt_volts = 10^((Pt_dBm - 30)/20);

% Thermal Noise
k = 1.38e-23;
T = 290;
B = (2.480 - 2.402)*1e9;
Pn = k*T*B;
Pn_dBW = 10*log10(Pn);
Pn_dBm = Pn_dBW + 30;
noise = wgn(1,L,Pn_dBW,'complex');

% Data Subcarriers
data_channel_freq = 1e6*[2402+2:2:2426-2  2426+2:2:2480-2]';
x_data = Pt_volts * exp(1j*2*pi.*data_channel_freq.*t) + noise;

% Advertisement Subcarriers
adv_channel_freq = 1e6*[2402 2426 2480]';
x_adv = Pt_volts * exp(1j*2*pi.*adv_channel_freq.*t) + noise;

% Phase Noise Specifications
phase_noise_freq = [ 1e3, 10e3, 100e3, 1e6, 10e6 ]; % Offset From Carrier
phase_noise_power = [ -84, -100, -96, -109, -122 ]; % Phase Noise power

% Add Phase Noise to Signal (Used function from internet)
x_adv = add_phase_noise( x_adv, Fs, phase_noise_freq, phase_noise_power );
x_data = add_phase_noise( x_data, Fs, phase_noise_freq, phase_noise_power );

% FFT for viewing channels in frequency domain
NFFT = length(x_data);
w = gausswin(NFFT);

x_data_tot = sum(x_data);
X_data = fftshift(fft(x_data_tot.*w'));
x_adv_tot = sum(x_adv);
X_adv = fftshift(fft(x_adv_tot.*w'));

% Normalize by FFT gain and add power from negative frequencies (*2)
X_data = 2*X_data/NFFT;                 
X_adv = 2*X_adv/NFFT;                 

f = Fs*(-0.5:1/NFFT:0.5-1/NFFT);

% Plot Channels
figure(1)
plot(f/1e9,20*log10(abs(X_data)) + 30,'k')
hold on
plot(f/1e9,20*log10(abs(X_adv)) + 30,'r')
hold off
xlim([2.4 2.482])
ylim([-60, Pt_dBm + 10])
xlabel('Frequency (GHz)')
ylabel('Power (dBm)')
title('Bluetooth Low Energy GFSK Channels')
legend('Data Subcarriers','Advertisment Subcarriers')

%% Get Power Generated from Piezo

[data1,data2] = arduinoData();              % This function gets the data collected from experiment

% Convert Data reading to Power Generated
expData_bitval = data2(900:end-400);        % Zoom in on region of interest
valRes = 1e6;                               % Resistor used to draw power in testing
valVoltage = 5;                             % Voltage held during test
adcBits = 2^10;                             % # ADC bits used
mapping = (valVoltage / adcBits);           % Mapping scalar to get bits to voltage
expData_Volt = expData_bitval * mapping;    % Scale to voltage
expData_Power = expData_Volt.^2 / valRes;   % Calculate Power Generated from V^2/R

% Create time vector
baud = 9600;                                % Baud rate of ADC
bytes_per_sec = baud / 8;                   % Convert to bytes/sec
t_piezo = (0:length(expData_bitval)-1) ./ bytes_per_sec * 2; % one step vs. two step

% Calculate power generated from one step (Averaged)
index_oneStep = 250;                        % Index where first step ended
avgPower_step = mean(expData_Power);
t_oneStep = t_piezo(1:index_oneStep);
duration_of_Step = t_oneStep(end);

% Plot Step Power Generated
figure(3)
plot(t_piezo,expData_Power*1e6,'k')
ylabel('Power Generated (uW)')
xlabel('time (s)')
title('Power Generated Over Time (Showing 5 Steps)')
grid on

%% Advertisment Packet Generation

inc = 1;
advertisement_Data_bytes = 1:31;
num_steps_needed = zeros(1,length(advertisement_Data_bytes));

% Loop over different packet lengths
for advertisement_Data_bytes_inc = advertisement_Data_bytes

    advDelay = 0.02;                        % Limited to 20ms to 10.24s
    n = 3;                                  % number of random numbers
    advRandom = 0 + rand(1,n) * (.01 - 0);  % random delay 0ms to 10ms
    advInterval = advDelay + advRandom;

    advAddress_bytes = 6;
    header_bytes = 2;

    payload_bytes = advAddress_bytes + advertisement_Data_bytes_inc;
    packetSize_bytes = header_bytes + payload_bytes; 
    packetSize_bits = packetSize_bytes * 8; %bytes x 8 bits
    bitRate = 1e6;                          %1 Mbps
    advPacket_length = (packetSize_bits) / bitRate;

    % Power Consumption for number of data bytes
    
    %See reference: https://www2.informatik.hu-berlin.de/~hs/Lehre/2014-WS_SP-CPS/20141028_Low-Power-Wireless-Standards_1.pdf
    %0.153 uWatts/bit for 0dBm Pt (taken from reference ^)
    power_per_bit = .153e-6;    
    power_per_packet = power_per_bit * packetSize_bits;

    energy_per_Packet_Joules = power_per_packet * advPacket_length;
    energy_per_Broadcast_Joules = 3 * energy_per_Packet_Joules; % Same packet 3 times on 3 adv. channels

    max_Rectifier_Efficiency = 40.6 / 100; %40.6 Reference --> https://www.brighthubengineering.com/consumer-appliances-electronics/96645-efficiency-of-ac-rectifiers/
    energy_per_Step_Joules = max_Rectifier_Efficiency * avgPower_step * duration_of_Step;

    % Ratio of Joules needed per broadcast to Joules generated per step
    num_steps_needed(inc) = energy_per_Broadcast_Joules/energy_per_Step_Joules;
    inc = inc + 1;

end

% Plot # of steps needed to send packet of given length
figure(4)
plot(advertisement_Data_bytes,num_steps_needed,'ro')
title('Number of Steps Needed to Power a Broadcast')
ylabel('# of Steps Needed')
xlabel('Packet Length (Bytes)')
grid on

% Adv Interval Required
walkingPace = (0.83:.2:1.95).';               % steps per second (Hz)
adv_Interval_Required = (2./walkingPace) * num_steps_needed; % Factor of two from only using one foot to generate power 

figure(5)
plot(advertisement_Data_bytes,adv_Interval_Required)
max_Interval_Allowed = 10.24;
hold on
plot([0 35],[max_Interval_Allowed max_Interval_Allowed],'r--')
title('Advertisement Interval Required vs. Packet Length')
ylabel('Advertisement Interval Required (s)')
xlabel('Packet Length (Bytes)')

for i = 1:length(walkingPace)
    str(i,:) = sprintf('%f',(walkingPace(i)));
end

leg = legend(str,'Location','NorthWest');
title(leg,'Walking Frequency (Hz)')
grid on
hold off
    
    

