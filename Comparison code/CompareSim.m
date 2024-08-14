%Compare data

trapSatData = allData(:,:,2);
interactionInt = allData(:,:,4);

%% low trap number
idx = [1 5 10];
yl = [0 35000];
xl = [0 1000];

figure
subplot(3,2,1)
plot(trapSatData(idx(1),:),'r')
ylim(yl)
xlim(xl)

subplot(3,2,2)
plot(interactionInt(idx(1),:),'b')
ylim(yl)
xlim(xl)

subplot(3,2,3)
plot(trapSatData(idx(2),:),'r')
ylim(yl)
xlim(xl)

subplot(3,2,4)
plot(interactionInt(idx(2),:),'b')
ylim(yl)
xlim(xl)

subplot(3,2,5)
plot(trapSatData(idx(3),:),'r')
ylim(yl)
xlim(xl)

subplot(3,2,6)
plot(interactionInt(idx(3),:),'b')
ylim(yl)
xlim(xl)
%% high trap number
idx = [9 13 17];

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
%% stats figure

figure
subplot(1,5,1)
plot(statsTS.nTraps,statsTS.med,'r')
hold on
plot(statsTS.nTraps,statsTS.medNN,'--r')
plot(statsTS.nTraps,statsInt.med,'b')
plot(statsTS.nTraps,statsInt.medNN,'--b')
xlabel('Trap concentration (trap/Vol)')
ylabel('Median signal (counts)')
axis square
title('Median signal')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,2)
plot(statsTS.nTraps,statsTS.width,'r')
hold on
plot(statsTS.nTraps,statsTS.widthNN,'--r')
plot(statsTS.nTraps,statsInt.width,'b')
plot(statsTS.nTraps,statsInt.widthNN,'--b')

axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Std')
title('Standard deviation')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,3)
plot(statsTS.nTraps,statsTS.width./statsTS.med,'r')
hold on
plot(statsTS.nTraps,statsTS.widthNN./statsTS.medNN,'--r')
plot(statsTS.nTraps,statsInt.width./statsInt.med,'b')
plot(statsTS.nTraps,statsInt.widthNN./statsInt.medNN,'--b')
axis square
ylim([0 0.5])
xlabel('Trap concentration (trap/Vol)')
ylabel('Std/med')
title('Standard deviation / median')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,4)
plot(statsTS.nTraps,statsTS.absAmp,'r')
hold on
plot(statsTS.nTraps,statsTS.absAmpNN,'--r')

plot(statsTS.nTraps,statsInt.absAmp,'b')
plot(statsTS.nTraps,statsInt.absAmpNN,'--b')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Absolute amplitude')
title('Abs. Amplitude')
legend({'no trap interaction','no noise','trap interaction','no noise'})

subplot(1,5,5)
plot(statsTS.nTraps,statsTS.relAmp,'r')
hold on
plot(statsTS.nTraps,statsTS.relAmpNN,'--r')
plot(statsTS.nTraps,statsInt.relAmp,'b')
plot(statsTS.nTraps,statsInt.relAmpNN,'--b')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Relative amplitude')
title('Relative amplitude')
legend({'no trap interaction','no noise','trap interaction','no noise'})
ylim([0 1])