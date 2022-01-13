

%% Input


%%



x = linspace(0,100*pi,10000);


%% Amplitude scan

freq = 10;

y = zeros(10, length(x));

amp = [0.0001 0.001 0.01 0.1 1 10 100 1000 10000 100000];
for i = 1:10
    
    y(i,:) = amp(i)*sin(freq*x);
    
    
    
end

%% Freq Scan

y = zeros(10, length(x));

freq = [0.01, 0.05, 0.1, 0.5, 1, 10,50];
amp = 1;
for i = 1:7
    
    y(i,:) = amp*sin(freq(i)*x);
    
    
    
end

