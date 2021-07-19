% Representation of two uncoded pulses
pulse1 = [1 1 1 1 1 1 1 1 1 1 1 1 1];
pulse2 = [1 1 1 1 1 1 1 1 1 1 1 1 1];

% Matched filter output (pulse correlated with itself)
correlator_output1 = xcorr(pulse1);
correlator_output2 = xcorr(pulse2);

% 13-baud Barker Code
c = 3e8; 
i_baud = -12:1:12;     % x-index, peak centered at 0
baud_length_t = 10e-6   % baud length in seconds
baud_length_km = c*baud_length_t/1e3; %baud length in km
i_baud_km=i_baud*baud_length_km;    % x-index in km

% Plot output power vs range for each pulse
figure(1); 
hold off;
plot(i_baud_km,correlator_output1); 
grid on;
xlabel('Range (km)');
ylabel('Output Power (arbitrary unites)');
hold on;
plot(i_baud_km,correlator_output2); grid on;
title('Matched Filter Output');