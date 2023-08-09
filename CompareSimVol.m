%Compare data

trapSatData = allData(:,:,2);
interactionInt = allData(:,:,4);

%% low trap number
idx = [1 3 4];

figure
subplot(3,2,1)
plot(trapSatData(idx(1),:),'r')
ylim([0 max(trapSatData(idx(1),:))])
xlim([0 1000])

subplot(3,2,2)
plot(interactionInt(idx(1),:),'b')
ylim([0 max(interactionInt(idx(1),:))])
xlim([0 1000])

subplot(3,2,3)
plot(trapSatData(idx(2),:),'r')
ylim([0 max(trapSatData(idx(2),:))])
xlim([0 1000])

subplot(3,2,4)
plot(interactionInt(idx(2),:),'b')
ylim([0 max(interactionInt(idx(2),:))])
xlim([0 1000])

subplot(3,2,5)
plot(trapSatData(idx(3),:),'r')
ylim([0 max(trapSatData(idx(3),:))])
xlim([0 1000])

subplot(3,2,6)
plot(interactionInt(idx(3),:),'b')
ylim([0 max(interactionInt(idx(3),:))])
xlim([0 1000])
%% high trap number
idx = [5 7 8];

figure
subplot(3,2,1)
plot(trapSatData(idx(1),:),'r')
ylim([0 max(trapSatData(idx(1),:))])
xlim([0 1000])

subplot(3,2,2)
plot(interactionInt(idx(1),:),'b')
ylim([0 max(interactionInt(idx(1),:))])
xlim([0 1000])

subplot(3,2,3)
plot(trapSatData(idx(2),:),'r')
ylim([0 max(trapSatData(idx(2),:))])
xlim([0 1000])

subplot(3,2,4)
plot(interactionInt(idx(2),:),'b')
ylim([0 max(interactionInt(idx(2),:))])
xlim([0 1000])

subplot(3,2,5)
plot(trapSatData(idx(3),:),'r')
ylim([0 max(trapSatData(idx(3),:))])
xlim([0 1000])

subplot(3,2,6)
plot(interactionInt(idx(3),:),'b')
ylim([0 max(interactionInt(idx(3),:))])
xlim([0 1000])
%% stats figure

figure
subplot(1,5,1)
hold on
plot(statsTS.nTraps,statsInt.med,'b')
plot(statsTS.nTraps,statsInt.medNN,'--b')
xlabel('Volume (a.u.)')
ylabel('Median signal (counts)')
axis square
title('Median signal')
legend({'trap interaction','no noise'})

subplot(1,5,2)
hold on
plot(statsTS.nTraps,statsInt.width,'b')
plot(statsTS.nTraps,statsInt.widthNN,'--b')

axis square
xlabel('Volume (a.u.)')
ylabel('Std')
title('Standard deviation')
legend({'trap interaction','no noise'})

subplot(1,5,3)
hold on
plot(statsTS.nTraps,statsInt.width./statsInt.med,'b')
plot(statsTS.nTraps,statsInt.widthNN./statsInt.medNN,'--b')
axis square
ylim([0 0.5])
xlabel('Volume (a.u.)')
ylabel('Std/med')
title('Standard deviation / median')
legend({'trap interaction','no noise'})

subplot(1,5,4)
hold on
plot(statsTS.nTraps,statsInt.absAmp,'b')
plot(statsTS.nTraps,statsInt.absAmpNN,'--b')
axis square
xlabel('Volume (a.u.)')
ylabel('Absolute amplitude')
title('Abs. Amplitude')
legend({'trap interaction','no noise'})

subplot(1,5,5)
hold on
plot(statsTS.nTraps,statsInt.relAmp,'b')
plot(statsTS.nTraps,statsInt.relAmpNN,'--b')
axis square
xlabel('Volume (a.u.)')
ylabel('Relative amplitude')
title('Relative amplitude')
legend({'trap interaction','no noise'})
ylim([0 1])