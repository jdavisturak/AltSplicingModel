
load('Analytic_Function_Theta5_2.mat');
close all
% Simple tests - heatmap time!!!

figure(1)
for f=1:5
    subplot(2,3,f)
    PARAMS = [1.0001 1.00015 f*1.0002/10,1/3, 1/30, 1/3];
    
    dat1  = zeros(100,100);
    
    % Default value:
    fprintf('Default value: %f \n', F_Theta5_2(PARAMS(1),PARAMS(2),PARAMS(3),PARAMS(4),PARAMS(5),PARAMS(6)))
    
    multsR = logspace(-3,3,size(dat1,1));
    multsS = logspace(-3,3,size(dat1,2));
    
    for ms = 1:length(multsS)
        for mr = 1:length(multsR)
            dat1(ms, mr) = F_Theta5_2(PARAMS(1)*multsR(mr),PARAMS(2)*multsS(ms),PARAMS(3),PARAMS(4),PARAMS(5),PARAMS(6));
        end
    end
    
    max(max(dat1))
    
    image(dat1,'CDataMapping','scaled'); set(gca,'CLim',[0 1]); colorbar
end

figure(2)
% Defaults
PARAMS = [1.0001 1.00015 1.0002/2,1/3, 1/30, 1/3];

dat2  = zeros(100,100);

% Default value:
fprintf('Default value: %f \n', F_Theta5_2(PARAMS(1),PARAMS(2),PARAMS(3),PARAMS(4),PARAMS(5),PARAMS(6)))

multsR = logspace(-1,1,size(dat2,1));
multsS = logspace(-1,1,size(dat2,2));

for ms = 1:length(multsS)
    for mr = 1:length(multsR)
        dat2(ms, mr) = F_Theta5_2(PARAMS(1)*multsR(mr),PARAMS(2)*multsR(mr),PARAMS(3)*multsS(ms),PARAMS(4),PARAMS(5),PARAMS(6));
    end
end

max(max(dat2))

image(dat2,'CDataMapping','scaled'); set(gca,'CLim',[0.5 1]); colorbar












