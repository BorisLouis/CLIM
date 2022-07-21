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
subplot(1,5,1)
plot(statsTS.med,'r')
hold on
plot(statsTS.medNN,'--r')
plot(statsInt.med,'b')
plot(statsInt.medNN,'--b')
xlabel('Trap concentration (trap/Vol)')
ylabel('Median signal (counts)')
axis square
title('Median signal')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,2)
plot(statsTS.width,'r')
hold on
plot(statsTS.widthNN,'--r')
plot(statsInt.width,'b')
plot(statsInt.widthNN,'--b')

axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Std')
title('Standard deviation')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,3)
plot(statsTS.width./statsTS.med,'r')
hold on
plot(statsTS.widthNN./statsTS.medNN,'--r')
plot(statsInt.width./statsInt.med,'b')
plot(statsInt.widthNN./statsInt.medNN,'--b')
axis square
ylim([0 0.5])
xlabel('Trap concentration (trap/Vol)')
ylabel('Std/med')
title('Standard deviation / median')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,4)
plot(statsTS.absAmp,'r')
hold on
plot(statsTS.absAmpNN,'--r')

plot(statsInt.absAmp,'b')
plot(statsInt.absAmpNN,'--b')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Absolute amplitude')
title('Abs. Amplitude')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,5)
plot(statsTS.relAmp,'r')
hold on
plot(statsTS.relAmpNN,'--r')
plot(statsInt.relAmp,'b')
plot(statsInt.relAmpNN,'--b')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Relative amplitude')
title('Relative amplitude')
legend({'no trap interaction','no noise','trap interaction','no noise'})
ylim([0 1])