%408 Assignment 3

%% Part 1 A
sample1 = importdata('408_ECG_3O.txt');
sample2 = importdata('408_ECG_3Y.txt');

%Canceling the 8 V offset in sample2
for i= 1:size(sample2)
     sample2(i,1) = sample2(i,1) - 8;
end

%Selecting 5s intervals from data points
for i = 1:1250
    interval1(i,1) = sample1(i,1);
    interval2(i,1) = sample2(i,1);
end

for i = (size(sample1,1)-1250):size(sample1,1)
    interval3(i,1) = sample1(i,1);
end

for i = (size(sample2,1)-1250):size(sample2,1)
    
    interval4(i,1) = sample2(i,1);
end

% % %Plot data points
% figure
% plot(interval1)
% title('408 ECG 3O First 5s')
% xlabel('Sample #')
% ylabel('Voltage')
% % 
% figure
% plot(interval2)
% title('408 ECG 3Y First 5s')
% xlabel('Sample #')
% ylabel('Voltage')
% 
% figure
% plot(interval3)
% title('408 ECG 3O Final 5s')
% xlabel('Sample #')
% ylabel('ECG Data')
% xlim([438260 439510])
% % 
% figure
% plot(interval4)
% title('408 ECG 3Y Final 5s')
% xlabel('Sample #')
% ylabel('Voltage')
% xlim([450846 452096])

%% Part 1 B :Amplitude Spectra
% %calculating ASD
% f = 250;
% figure
% tiledlayout(2,1)
% nexttile
%     psdest = psd(spectrum.periodogram,interval1,'Fs',f);
%     plot(psdest.Frequencies,psdest.Data);
%     xlabel('Hz'); grid on;
%       xdft1 = fft(interval1);
%   % assume x is even length
%   xdft1 = xdft1(1:length(interval1)/2+1);
%   freq = 0:f/length(interval1):f/2;
%   plot(freq,abs(xdft1));
%   title('Amplitude Spectral Density 408 ECG 3O')
%   xlabel('Hz');
% nexttile
%     psdest = psd(spectrum.periodogram,interval2,'Fs',f);
%     plot(psdest.Frequencies,psdest.Data);
%     xlabel('Hz'); grid on;
%       xdft2 = fft(interval2);
%   % assume x is even length
%   xdft2 = xdft2(1:length(interval2)/2+1);
%   freq = 0:f/length(interval2):f/2;
%   plot(freq,abs(xdft2));
%   title('Amplitude Spectral Density 408 ECG 3Y')
%   xlabel('Hz');

%% Part 1 C: Intervals between heart beats

timing1 = zeros(439510,1);
timing2 = zeros(452096,1);

for i= 1:size(sample1(:,1))
    if sample1(i,1) > 1.5
        heartbeat1(i,1) = sample1(i,1);
        timing1(i,1) = i;
    else
        heartbeat1(i,1) = 0;
    end
end

for i= 1:size(sample2(:,1))
    if sample2(i,1) > 1.8
        heartbeat2(i,1) = sample2(i,1);
        timing2(i,1) = i;
 
    else
        heartbeat2(i,1) = 0;
    end
end

%Get NN interval
s1 = [heartbeat1(:,1) timing1(:,1)];
s1_1 = circshift(s1,1);
s2 = [heartbeat2(:,1) timing2(:,1)];
s2_1 = circshift(s2,1);
beat1 = s1-s1_1;
beat2 = s2-s2_1;
j = 1;
t = 1;

for i= 1:size(sample1(:,1))
    if beat1(i,2) > 1
        beat_1(j,1) = beat1(i,2);
        j = j + 1;
    end
end
for i= 1:size(sample2(:,1))
    if beat2(i,2) > 1
        beat_2(t,1) = beat2(i,2);
        t = t + 1;
    end
end

beat2_3 = circshift(beat_1,1);
beat2_4 = circshift(beat_2,1);

beat1_diff = beat_1 - beat2_3;
beat2_diff = beat_2 - beat2_4;

%Get first 5 min of both samples
for i = 1:size(beat_1)
    if (beat_1(i) <= 75000)
        b_1(i) = beat_1(i);
    end
end
for i = 1:size(beat_2)
    if (beat_2(i) <= 75000)
        b_2(i) = beat_2(i);
    end
end
beat1_diff_5 = b_1 - circshift(b_1,1);
beat1_diff_5(beat1_diff_5 < 0) = [];
beat2_diff_5 = b_2 - circshift(b_2,1);
beat2_diff_5(beat2_diff_5 < 0) = [];

beat1_diff(beat1_diff < 0) = []; %Delete the one off data point from above
beat2_diff(beat2_diff < 0) = [];

%% Part 1 D Tachograms
% figure
% plot(beat1_diff_5)
% title('Tachogram for 408 ECG 3O')
% xlabel('Beat #')
% ylabel('Time between R Waves')
% % xlim([0 75000])
% % 
% figure
% plot(beat2_diff_5)
% title('Tachogram for 408 ECG 3Y')
% xlabel('Beat #')
% ylabel('Time between R Waves')
% % xlim([0 75000])

%% Part 1 E Mean and Range
count = 1;

for i=1:round(1492/6):1492-round(1492/6) %Sample 1
meanNN1(count)= mean(beat1_diff(i:i+(1492/6)));
rangeNN1(count)= range(beat1_diff(i:i+(1492/6)));
count=count+1;
end
count = 1;
for i=1:round(1706/6):1706-round(1706/6) %Sample 2
meanNN2(count)= mean(beat2_diff(i:i+(1706/6)));
rangeNN2(count)= range(beat2_diff(i:i+(1706/6)));
count=count+1;
end
% 
% % %Plot mean and range comparison for sample 1
% figure
% plot(meanNN1)
% title('NN Interval Mean for 408 ECG 3O')
% xlabel('Segment #')
% ylabel('NN Interval Time')
% 
% figure
% plot(rangeNN1)
% title('NN Interval Range for 408 ECG 3O')
% xlabel('Segment #')
%  ylabel('NN Interval Time')
% 
% %Plot mean and range comparison for sample 2
% figure
% plot(meanNN2)
% title('NN Interval Mean for 408 ECG 3Y')
% xlabel('Segment #')
%  ylabel('NN Interval Time')
% 
% figure
% plot(rangeNN2)
% title('NN Interval Range for 408 ECG 3Y')
% xlabel('Segment #')
%  ylabel('NN Interval Time')

%% Part 2 

%Calculate the standard deviation
count = 1;

for i=1:round(1492/6):1492-round(1492/6) %Sample 1
stdNN1(count)= std(beat1_diff(i:i+(1492/6)));
count=count+1;
end
count = 1;
for i=1:round(1706/6):1706-round(1706/6) %Sample 2
stdNN2(count)= std(beat2_diff(i:i+(1706/6)));
count=count+1;
end
% 
% figure
% plot(stdNN1)
% title('NN Interval STD for Sample 1')
% xlabel('Segment #')
% ylabel('STD')
% 
% figure
% plot(stdNN2)
% title('NN Interval STD for Sample 2')
% xlabel('Segment #')
% ylabel('STD')

%% Part 3
power_spect1 = 1/(439510.^ 2) * abs(fft(sample1).^2); %https://www.dsprelated.com/showarticle/1004.php
% power_spect2 = 1/(452096.^ 2) * (xdft2.^2);
% 
% %Plot the Power Spectrums
% figure
% plot(power_spect1)
% title('Power Spectrum 1')
% xlabel('Time')
% ylabel('Power')
% 
% figure
% plot(power_spect2)
% title('Power Spectrum 2')
% xlabel('Time')
% ylabel('Power')

%% Part 4
%Poincare plot
beat1_diff_1 = circshift(beat1_diff,1);

figure
plot(beat1_diff, beat1_diff_1)
title('Poincare Plot 408 ECG 3O')
xlabel('NN(n)')
ylabel('NN(n+1)')

beat2_diff_1 = circshift(beat2_diff,1);

figure
plot(beat2_diff, beat2_diff_1)
title('Poincare Plot 408 ECG 3Y')
xlabel('NN(n)')
ylabel('NN(n+1)')

std1 = std(beat1_diff)
std2 = std(beat2_diff)

ratio_std = std1/std2