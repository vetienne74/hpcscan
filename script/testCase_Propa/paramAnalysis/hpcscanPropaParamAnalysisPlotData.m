
close all ; clear all ;

DIR  = '.' ;
% log file name with .log extension
%FILE = 'hpcscanPropaParamAnalysisShaheen' ;
FILE = 'hpcscan.perf.Propa' ;

pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

valTime     = val.data(:,14) ;
valN1       = val.data(:,6) ;
valError    = val.data(:,15) ;
valOrder    = val.data(:,9) ;
valGpoint   = val.data(:,12) ;
valGbyte    = val.data(:,13) ;

valDt       = val.data(:,17) ;
valStableDt = val.data(:,18) ;

sizeVal = size(val.data) ;
nConfig = sizeVal(1)

figure('Position',[100 100 1100 900])

% plot Error versus Time
subplot(2,2,2)
hold on; grid on;
xlabel('Elapse Time (s)')
ylabel('L1 Error')
title('L1 Error vs Elapse Time', 'FontSize', 12)

maxAllowedError = 0.01 ;
elapseTimeOptimal = 9.99e+9;

for ii=1:nConfig
    
    % search for optimal config
    if valError(ii) < maxAllowedError
        if valTime(ii) < elapseTimeOptimal
            iConfigOptimal = ii ;
            elapseTimeOptimal = valTime(ii) ;
        end
    end
    
    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        if valOrder(ii) == 2
            plot(valTime(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valTime(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valTime(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'm', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valTime(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2)
        else
            plot(valTime(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'LineWidth', 2)
        end
    else
        if valOrder(ii) == 2
            plot(valTime(ii), valError(ii), '+', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valTime(ii), valError(ii), '+', 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valTime(ii), valError(ii), '+', 'MarkerEdgeColor', 'm', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valTime(ii), valError(ii), '+', 'MarkerEdgeColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        else
            plot(valN1(ii), valError(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% represent optimal config with a star
plot(valTime(iConfigOptimal), valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

% display info on best best config
xText = 0.7*min(valTime(:)) ;
yText = 1.4*valError(iConfigOptimal) ;
yTextDelta = yText*0.36 ;
text(xText, yText, 0, sprintf('Target Err. %2.2f',maxAllowedError))
text(xText, yText-yTextDelta, 0, 'Optim. config.')
text(xText, yText-1.65*yTextDelta, 0, sprintf('FD O%i %s',valOrder(iConfigOptimal),...
    string(val.textdata(iConfigOptimal,4))))
ratioCFL= valDt(iConfigOptimal) / valStableDt(iConfigOptimal) * 100;
text(xText, yText-2.05*yTextDelta, 0, sprintf('N=%i, %2.2f%% CFL',valN1(iConfigOptimal), ratioCFL))

ax=gca
ax.XScale='log'
ax.YScale='log'

% plot Error versus N1
subplot(2,2,1)
hold on; grid on;
xlabel('N')
ylabel('L1 Error')
title('L1 Error vs N', 'FontSize', 12)
for ii=1:nConfig
    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        if valOrder(ii) == 2
            plot(valN1(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN1(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN1(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'm', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN1(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2)
        else
            plot(valN1(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'LineWidth', 2)
        end
    else
        if valOrder(ii) == 2
            plot(valN1(ii), valError(ii), '+', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN1(ii), valError(ii), '+', 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN1(ii), valError(ii), '+', 'MarkerEdgeColor', 'm', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN1(ii), valError(ii), '+', 'MarkerEdgeColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        else
            plot(valN1(ii), valError(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% represent optimal config with a star
plot(valN1(iConfigOptimal), valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
ax.YScale='log'

% plot Gpoint versus N1
subplot(2,2,3)
hold on; grid on;
xlabel('N')
ylabel('GPoint/s')
title('GPoint/s vs N', 'FontSize', 12)
for ii=1:nConfig
    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        if valOrder(ii) == 2
            plot(valN1(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN1(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN1(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'm', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN1(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2)
        else
            plot(valN1(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'LineWidth', 2)
        end
    else
        if valOrder(ii) == 2
            plot(valN1(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN1(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN1(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'm', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN1(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        else
            plot(valN1(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% represent optimal config with a star
plot(valN1(iConfigOptimal), valGpoint(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
ax.YScale='log'

% plot Gbyte versus N1
subplot(2,2,4)
hold on; grid on;
xlabel('N')
ylabel('GByte/s')
title('GByte/s vs N', 'FontSize', 12)
for ii=1:nConfig
    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        if valOrder(ii) == 2
            plot(valN1(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN1(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN1(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'm', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN1(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2)
        else
            plot(valN1(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'LineWidth', 2)
        end
    else
        if valOrder(ii) == 2
            plot(valN1(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN1(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN1(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'm', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN1(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        else
            plot(valN1(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% represent optimal config with a star
plot(valN1(iConfigOptimal), valGbyte(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
ax.YScale='log'

% save figure
figName = sprintf('%s.jpg', FILE) ;
print(figName, '-djpeg')


% END
