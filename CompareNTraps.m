

allData = [T500.statsInt;T1000.statsInt;T2000.statsInt];

allData.nTraps = [500; 1000; 2000];


figure
subplot(1,5,1)
plot(allData.nTraps, allData.med,'r')
hold on
plot(allData.nTraps, allData.medNN,'--r')
xlabel('Trap concentration (trap/Vol)')
ylabel('Median signal (counts)')
axis square
title('Median signal')
legend({'trap interaction','no noise'})

subplot(1,5,2)
plot(allData.nTraps, allData.width,'r')
hold on
plot(allData.nTraps, allData.widthNN,'--r')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Std')
title('Standard deviation')
legend({'trap interaction','no noise'})

subplot(1,5,3)
plot(allData.nTraps, allData.width./allData.med,'r')
hold on
plot(allData.nTraps, allData.widthNN./allData.medNN,'--r')
axis square
ylim([0 0.5])
xlabel('Trap concentration (trap/Vol)')
ylabel('Std/med')
title('Standard deviation / median')
legend({'trap interaction','no noise'})

subplot(1,5,4)
plot(allData.nTraps, allData.absAmp,'r')
hold on
plot(allData.nTraps, allData.absAmpNN,'--r')

axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Absolute amplitude')
title('Abs. Amplitude')
legend({'trap interaction','no noise'})

subplot(1,5,5)
plot(allData.nTraps, allData.relAmp,'r')
hold on
plot(allData.nTraps, allData.relAmpNN,'--r')
axis square
xlabel('Trap concentration (trap/Vol)')
ylabel('Relative amplitude')
title('Relative amplitude')
legend({'trap interaction','no noise'})
ylim([0 1])