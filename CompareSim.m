%Compare data

trapSatData = allData(:,:,1);
interactionInt = allData(:,:,2);

%% low trap number
idx = [1 5 10];

figure
subplot(3,2,1)
plot(trapSatData(idx(1),:),'r')
ylim([0 12000])
xlim([0 1000])

subplot(3,2,2)
plot(interactionInt(idx(1),:),'b')
ylim([0 12000])
xlim([0 1000])

subplot(3,2,3)
plot(trapSatData(idx(2),:),'r')
ylim([0 12000])
xlim([0 1000])

subplot(3,2,4)
plot(interactionInt(idx(2),:),'b')
ylim([0 12000])
xlim([0 1000])

subplot(3,2,5)
plot(trapSatData(idx(3),:),'r')
ylim([0 12000])
xlim([0 1000])

subplot(3,2,6)
plot(interactionInt(idx(3),:),'b')
ylim([0 12000])
xlim([0 1000])
%% high trap number
idx = [10 20 50];

figure
subplot(3,2,1)
plot(trapSatData(idx(1),:),'r')
%ylim([0 12000])
xlim([0 1000])

subplot(3,2,2)
plot(interactionInt(idx(1),:),'b')
%ylim([0 12000])
xlim([0 1000])

subplot(3,2,3)
plot(trapSatData(idx(2),:),'r')
%ylim([0 12000])
xlim([0 1000])

subplot(3,2,4)
plot(interactionInt(idx(2),:),'b')
%ylim([0 12000])
xlim([0 1000])

subplot(3,2,5)
plot(trapSatData(idx(3),:),'r')
%ylim([0 12000])
xlim([0 1000])

subplot(3,2,6)
plot(interactionInt(idx(3),:),'b')
%ylim([0 12000])
xlim([0 1000])
%% stats figure

figure
subplot(1,4,1)
plot(statsTS.med,'r')
hold on
plot(statsTS.medNN,'--r')
plot(statsInt.med,'b')
xlabel('Trap concentration (trap/Vol)')
ylabel('Median signal (counts)')
axis square
title('Median signal')
legend({'no trap interaction','trap interaction'})

subplot(1,4,2)
plot(statsTS.width./statsTS.med,'r')
hold on
plot(statsInt.width./statsInt.med,'b')
axis square
ylim([0 0.5])
xlabel('Trap concentration (trap/Vol)')
ylabel('Std/med')
title('Standard deviation / median')
legend({'no trap interaction','trap interaction'})

subplot(1,4,3)
plot(statsTS.absAmp,'r')
hold on
plot(statsInt.absAmp,'b')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Absolute amplitude')
title('Abs. Amplitude')
legend({'no trap interaction','trap interaction'})

subplot(1,4,4)
plot(statsTS.relAmp,'r')
hold on
plot(statsInt.relAmp,'b')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Relative amplitude')
title('Relative amplitude')
legend({'no trap interaction','trap interaction'})

