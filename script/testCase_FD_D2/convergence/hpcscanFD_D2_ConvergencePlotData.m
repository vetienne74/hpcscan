
close all ; clear all ;

DIR  = '.' ;

% derivative order (1 or 2)
ORD = 2 ;

% log file name without .log extension
FILE = sprintf('hpcscan.perf.FD_D%i', ORD) ;

% derivative (1=x1, 2=x2, 3=x3, 4=sum/Laplacian)
DER = 4 ;

% define target error
maxAllowedError = 0.01 ;

% define hardware memory bandwith
memBwdth = 44.0 ;

% min and max for plot
minErrPlot = 1.0e-8 ;
maxErrPlot = 1.0 ;
minGPoint  = 0.0 ;
maxGPoint  = 5.0 ;
minGByte   = 0.0 ;
maxGByte   = 200.0 ;
minTime    = 0.001 ;
maxTime    = 10.0 ;


%pathFile = sprintf('%s/%s.log.fp64', DIR, FILE) ;
pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

valOrder    = val.data(:,9) ;

if (DER == 1)
    valN        = val.data(:,6) ;
    valGpoint   = val.data(:,12) ;
    valGbyte    = val.data(:,13) ;
    valTime     = val.data(:,14) ;
    valError    = val.data(:,15) ;
    DER_AXIS = 'Axis1' ;
elseif (DER == 2)
    valN        = val.data(:,7) ;
    valGpoint   = val.data(:,18) ;
    valGbyte    = val.data(:,19) ;
    valTime     = val.data(:,20) ;
    valError    = val.data(:,21) ;
    DER_AXIS = 'Axis2' ;
elseif (DER == 3)
    valN        = val.data(:,8) ;
    valGpoint   = val.data(:,24) ;
    valGbyte    = val.data(:,25) ;
    valTime     = val.data(:,26) ;
    valError    = val.data(:,27) ;
    DER_AXIS = 'Axis3' ;
elseif (DER == 4)
    valN        = val.data(:,8) ;
    valGpoint   = val.data(:,30) ;
    valGbyte    = val.data(:,31) ;
    valTime     = val.data(:,32) ;
    valError    = val.data(:,33) ;
    DER_AXIS = 'Laplacian' ;
end

sizeVal = size(val.data) ;
nConfig = sizeVal(1)


figure('Position',[100 100 1100 900])

maxErrX(1:nConfig) = val.data(:,6) ;
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
    
    plot(valTime(ii), valError(ii), 'sw', 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', 5, 'LineWidth', 2)     

end

% represent optimal config with a star
plot(valTime(iConfigOptimal), valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

% display info on best best config
TITLE = sprintf('L1 Error vs Elapse Time \n Optimal FD O%i with N=%i', ...
    valOrder(iConfigOptimal),valN(iConfigOptimal)) ;
title(TITLE, 'FontSize', 12)

ax=gca
ax.XScale='log'
ax.YScale='log'

ylim([minErrPlot maxErrPlot])
xlim([minTime maxTime])

%----------------------
% plot Error versus N1
%----------------------

subplot(2,2,1)
hold on; grid on;
xlabel('N')
ylabel('L1 Error')
title('L1 Error vs N', 'FontSize', 12)

TITLE = sprintf('L1 Error vs N \n Derivative D%i - %s ', ORD, DER_AXIS)
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
    
    plot(valN(ii), valError(ii), 'sw', 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', 5, 'LineWidth', 2)   

end

% plot horizontal line with allowed error
plot(maxErrX, maxErrY, '-k', 'LineWidth', 1.5)

% represent optimal config with a star
plot(valN(iConfigOptimal), valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
ax.YScale='log'

ylim([minErrPlot maxErrPlot])

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
    
    plot(valN(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', 5, 'LineWidth', 2) 

end

% represent optimal config with a star
plot(valN(iConfigOptimal), valGpoint(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
%ax.YScale='log'
ylim([minGPoint maxGPoint])

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
    
    plot(valN(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', [colorR colorG colorB] , ...
        'MarkerFaceColor', [colorR colorG colorB] , 'MarkerSize', 5, 'LineWidth', 2)  

end

% represent optimal config with a star
plot(valN(iConfigOptimal), valGbyte(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
%ax.YScale='log'

ylim([minGByte maxGByte])

memBwdthX(1:nConfig) = val.data(:,6) ;
memBwdthY(1:nConfig) = memBwdth ;

% plot horizontal line with hardware memory bandwdith
plot(memBwdthX, memBwdthY, '-k', 'LineWidth', 1.5)

% save figure
figName = sprintf('%s-%s.jpg', FILE, DER_AXIS) ;
print(figName, '-djpeg')

fprintf('Error min %g - max %g\n', min(valError(:)), max(valError(:)))
fprintf('Time min %g - max %g\n',  min(valTime(:)), max(valTime(:)))

% END
