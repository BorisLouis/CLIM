clc
close all
%% user input
P = {'OD2','OD2.5','OD3','OD3.3'};

%%

figure(1)
hold on
if exist('corrOutputP1')
    scatter([corrOutputP1.results.tau]/(2*pi),[corrOutputP1.results.beta],10,'filled','MarkerFaceAlpha',0.5);
end
if exist('corrOutputP2')
    scatter([corrOutputP2.results.tau]/(2*pi),[corrOutputP2.results.beta],10,'filled','MarkerFaceAlpha',0.5);
end
if exist('corrOutputP3')
    scatter([corrOutputP3.results.tau]/(2*pi),[corrOutputP3.results.beta],10,'filled','MarkerFaceAlpha',0.5);
end
if exist('corrOutputP4')
    scatter([corrOutputP4.results.tau]/(2*pi),[corrOutputP4.results.beta],10,'filled','MarkerFaceAlpha',0.5);
end


xlim([0.3 2])
ylim([1.5 2.5])
set(gca,'XScale','log')
axis square
box on

legend(P)
xlabel('Tau(s)')
ylabel('Beta')

%%
figure(2)
edges = [0:0.01:0.5];
hold on
if exist('corrOutputP1')
    histogram([corrOutputP1.results.Amp],edges);
end
if exist('corrOutputP2')
    histogram([corrOutputP2.results.Amp],edges);
end
if exist('corrOutputP3')
    histogram([corrOutputP3.results.Amp],edges);
end
if exist('corrOutputP4')
    histogram([corrOutputP4.results.Amp],edges);
end
axis square
box on
legend(P)
xlabel('Blinking Amplitude')
ylabel('Occurrence')

%% 

figure(3)
edges = [0:0.01:0.5];
hold on
if exist('corrOutputP1')
    histogram([corrOutputP1.results.meanSil],edges);
end
if exist('corrOutputP2')
    histogram([corrOutputP2.results.meanSil],edges);
end
if exist('corrOutputP3')
    histogram([corrOutputP3.results.meanSil],edges);
end
if exist('corrOutputP4')
    histogram([corrOutputP4.results.meanSil],edges);
end
axis square
box on
legend(P)
xlabel('SIlhouette')
ylabel('Occurrence')

