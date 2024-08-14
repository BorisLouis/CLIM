clc
clear
close all
%% Compare size

figure

subplot(1,3,1)
hold on
scatter([corrOutputBig.results.tau],[corrOutputBig.results.beta],10,'filled','MarkerFaceAlpha',0.5)
scatter([corrOutputMed.results.tau],[corrOutputMed.results.beta],10,'filled','MarkerFaceAlpha',0.5)
scatter([corrOutputSmall.results.tau],[corrOutputSmall.results.beta],10,'filled','MarkerFaceAlpha',0.5)
scatter([corrOutputPorous.results.tau],[corrOutputPorous.results.beta],10,'filled','MarkerFaceAlpha',0.5)
xlim([0.1 7])
ylim([1.5 2.5])
axis square
box on
title('Timescale')
legend({'L','M','S','P'})
xlabel('Tau (s)')
ylabel('Beta')

subplot(1,3,2)
hold on
histogram([corrOutputBig.results.Amp],20)
histogram([corrOutputMed.results.Amp],20)
histogram([corrOutputSmall.results.Amp],20)
histogram([corrOutputPorous.results.Amp],20)
axis square
box on
title('Amplitude')
xlabel('Relative Amplitude')
ylabel('Occurence')
legend({'L','M','S','P'})

subplot(1,3,3)
hold on
histogram([corrOutputBig.results.meanSil],50)
histogram([corrOutputMed.results.meanSil],50)
histogram([corrOutputSmall.results.meanSil],50)
histogram([corrOutputPorous.results.meanSil],50)
axis square
box on
title('Silhouette')
legend({'L','M','S','P'})
xlabel('Silhouette Score')
ylabel('Occurence')