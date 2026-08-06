
close all ; clear all ;

% define target error
maxAllowedError = 0.01 ;

% set axis limits for graphs
PLOT_MINMAX = 0 ; % 0=no minmax, 1=with minmax

if (PLOT_MINMAX == 1)
    minErrPlot = 1.0e-8 ;
    maxErrPlot = 1.0 ;
    minGPoint  = 0.0 ;
    maxGPoint  = 5.0 ;
    minGByte   = 0.0 ;
    maxGByte   = 200.0 ;
    minTime    = 0.001 ;
    maxTime    = 10.0 ;
    memBwdth = 44.0 ; % define hardware memory bandwith
end

% END OF USER INPUT PARAMETERS

% number of modes defined in test case
NMODE = 150 ;

% read data from log file
%------------------------
DIR  = '.' ;
FILE = 'hpcscan.perf.Propa' ;
pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

valTime     = val.data(:,14) ;
valN        = val.data(:,6) ;
valError    = val.data(:,15) ;
valOrder    = val.data(:,9) ;
valGpoint   = val.data(:,12) ;
valGbyte    = val.data(:,13) ;

valDt       = val.data(:,17) ;
valStableDt = val.data(:,18) ;

sizeVal = size(val.data) ;
nConfig = sizeVal(1)

figure('Position',[100 100 1100 900])

maxErrX(1:nConfig) = valN(:) ;
maxErrY(1:nConfig) = maxAllowedError ;

%------------------------
% plot Error versus Time
%------------------------

subplot(2,2,2)
hold on; grid on;
xlabel('Elapse Time (s)')
ylabel('L1 Error')

elapseTimeOptimal = 9.99e+9;

for ii=1:nConfig

    % search for optimal config
    if valError(ii) < maxAllowedError
        if valTime(ii) < elapseTimeOptimal
            iConfigOptimal = ii ;
            elapseTimeOptimal = valTime(ii) ;
        end
    end

    if (valOrder(ii) <= 8)
        colorR = 0 ;
        colorG = (valOrder(ii) - 1) / 7 ;
        colorB = 1 - colorG;
    else
        colorR = (valOrder(ii) - 7) / 9 ;
        colorG = 1 - colorR ;
        colorB = 0 ;
    end

    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        markerType = 'sw' ;
        markerSize = 5 ;
    else
        markerType = '+' ;
        markerSize = 10 ;
    end
    
    plot(valTime(ii), valError(ii), markerType, 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', markerSize, 'LineWidth', 2)     

end

% represent optimal config with a star
plot(valTime(iConfigOptimal), valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

% display info on best best config
ratioCFL= valDt(iConfigOptimal) / valStableDt(iConfigOptimal) * 100;
TITLE = sprintf('L1 Error vs Elapse Time \n Optimal %s O%i / N=%i / %2.2f%% CFL', ...
    string(val.textdata(iConfigOptimal,4)), valOrder(iConfigOptimal),valN(iConfigOptimal), ratioCFL) ;
title(TITLE, 'FontSize', 12)

ax=gca
ax.XScale='log'
ax.YScale='log'

if (PLOT_MINMAX == 1)
    ylim([minErrPlot maxErrPlot])
    xlim([minTime maxTime])
end

%---------------------------------------
% plot Error versus npoint / wavelength
%---------------------------------------

nlambda = NMODE / 2 ;

subplot(2,2,1)
hold on; grid on;
xlabel('# points / wavelength')
ylabel('L1 Error')
title('L1 Error vs N', 'FontSize', 12)

TITLE = 'L1 Error vs spatial sampling' ;
title(TITLE, 'FontSize', 12)

for ii=1:nConfig      

    if (valOrder(ii) <= 8)
        colorR = 0 ;
        colorG = (valOrder(ii) - 1) / 7 ;
        colorB = 1 - colorG;
    else
        colorR = (valOrder(ii) - 7) / 9 ;
        colorG = 1 - colorR ;
        colorB = 0 ;
    end

    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        markerType = 'sw' ;
        markerSize = 5 ;
    else
        markerType = '+' ;
        markerSize = 10 ;
    end
    
    plot(valN(ii)/nlambda, valError(ii), markerType, 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', markerSize, 'LineWidth', 2)   

end

% plot horizontal line with allowed error
plot(maxErrX/nlambda, maxErrY, '-k', 'LineWidth', 1.5)

% represent optimal config with a star
plot(valN(iConfigOptimal)/nlambda, valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
ax.YScale='log'

if (PLOT_MINMAX == 1)
    ylim([minErrPlot maxErrPlot])
end

%-----------------------
% plot Gpoint versus N1
%-----------------------

subplot(2,2,3)
hold on; grid on;
xlabel('N')
ylabel('GPoint/s')
title('GPoint/s vs N', 'FontSize', 12)
for ii=1:nConfig

    if (valOrder(ii) <= 8)
        colorR = 0 ;
        colorG = (valOrder(ii) - 1) / 7 ;
        colorB = 1 - colorG;
    else
        colorR = (valOrder(ii) - 7) / 9 ;
        colorG = 1 - colorR ;
        colorB = 0 ;
    end

    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        markerType = 'sw' ;
        markerSize = 5 ;
    else
        markerType = '+' ;
        markerSize = 10 ;
    end
    
    plot(valN(ii), valGpoint(ii), markerType, 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', markerSize, 'LineWidth', 2) 

end

% represent optimal config with a star
plot(valN(iConfigOptimal), valGpoint(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
%ax.YScale='log'

if (PLOT_MINMAX == 1)
    ylim([minGPoint maxGPoint])
end

%----------------------
% plot Gbyte versus N1
%----------------------

subplot(2,2,4)
hold on; grid on;
xlabel('N')
ylabel('GByte/s')
title('GByte/s vs N', 'FontSize', 12)
for ii=1:nConfig

    if (valOrder(ii) <= 8)
        colorR = 0 ;
        colorG = (valOrder(ii) - 1) / 7 ;
        colorB = 1 - colorG;
    else
        colorR = (valOrder(ii) - 7) / 9 ;
        colorG = 1 - colorR ;
        colorB = 0 ;
    end

    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        markerType = 'sw' ;
        markerSize = 5 ;
    else
        markerType = '+' ;
        markerSize = 10 ;
    end
    
    plot(valN(ii), valGbyte(ii), markerType, 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', markerSize, 'LineWidth', 2)  

end

% represent optimal config with a star
plot(valN(iConfigOptimal), valGbyte(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
%ax.YScale='log'

if (PLOT_MINMAX == 1)
    ylim([minGByte maxGByte])

    memBwdthX(1:nConfig) = val.data(:,6) ;
    memBwdthY(1:nConfig) = memBwdth ;

    % plot horizontal line with hardware memory bandwdith
    plot(memBwdthX, memBwdthY, '-k', 'LineWidth', 1.5)
end

% save figure
figName = sprintf('%s-4fig.jpg', FILE) ;
print(figName, '-djpeg')

%---------------------------------------
% Same but only one plot
% plot Error versus npoint / wavelength
%---------------------------------------

figure
hold on; grid on;
xlabel('# points / wavelength')
ylabel('L1 Error')
title('L1 Error vs N', 'FontSize', 12)

TITLE = 'L1 Error vs spatial sampling' ;
title(TITLE, 'FontSize', 12)

for ii=1:nConfig      

    if (valOrder(ii) <= 8)
        colorR = 0 ;
        colorG = (valOrder(ii) - 1) / 7 ;
        colorB = 1 - colorG;
    else
        colorR = (valOrder(ii) - 7) / 9 ;
        colorG = 1 - colorR ;
        colorB = 0 ;
    end

    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        markerType = 'sw' ;
        markerSize = 5 ;
    else
        markerType = '+' ;
        markerSize = 10 ;
    end
    
    plot(valN(ii)/nlambda, valError(ii), markerType, 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', markerSize, 'LineWidth', 2)   

end

% plot horizontal line with allowed error
plot(maxErrX/nlambda, maxErrY, '-k', 'LineWidth', 1.5)

% represent optimal config with a star
%plot(valN(iConfigOptimal)/nlambda, valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
ax.YScale='log'

if (PLOT_MINMAX == 1)
    ylim([minErrPlot maxErrPlot])
end

% save figure
figName = sprintf('%s.jpg', FILE) ;
print(figName, '-djpeg')

fprintf('Error min %g - max %g\n', min(valError(:)), max(valError(:)))
fprintf('Time min %g - max %g\n',  min(valTime(:)), max(valTime(:)))

% END
