[data1,data2] = arduinoData();

data2_volts = data2 * (5 / 1024);
data2_Power = data2_volts.^2 / 1e6;
I = data2_volts / 1e6;

current_supplied = data2_Power / 5;
figure(1)
plot(current_supplied*1e6)
ylabel('uA')