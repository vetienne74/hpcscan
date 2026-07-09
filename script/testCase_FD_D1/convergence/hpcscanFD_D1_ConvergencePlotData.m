
%
close all ; clear all ;

DIR  = '.' ;

% derivative order (1 or 2)
ORD = 1 ;

% log file name without .log extension
FILE = sprintf('hpcscan.perf.FD_D%i', ORD) ;

% derivative (1=x1, 2=x2, 3=x3)
DER = 3 ;

% define max allowed error
maxAllowedError = 0.01 ;

pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

valOrder    = val.data(:,9) ;

if (DER == 1)
    valN        = val.data(:,6) ;
    valGpoint   = val.data(:,12) ;
    valGbyte    = val.data(:,13) ;
    valTime     = val.data(:,14) ;
    valError    = val.data(:,15) ;
elseif (DER == 2)
    valN        = val.data(:,7) ;
    valGpoint   = val.data(:,18) ;
    valGbyte    = val.data(:,19) ;
    valTime     = val.data(:,20) ;
    valError    = val.data(:,21) ;
elseif (DER == 3)
    valN        = val.data(:,8) ;
    valGpoint   = val.data(:,24) ;
    valGbyte    = val.data(:,25) ;
    valTime     = val.data(:,26) ;
    valError    = val.data(:,27) ;
end



sizeVal = size(val.data) ;
nConfig = sizeVal(1)

figure('Position',[100 100 1100 900])

maxErrX(1:nConfig) = val.data(:,6) ;
maxErrY(1:nConfig) = maxAllowedError ;

% plot Error versus Time
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
            plot(valN(ii), valError(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% represent optimal config with a star
plot(valTime(iConfigOptimal), valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)


% display info on best best config

%xText = 0.7*min(valTime(:)) ;
xText = valTime(iConfigOptimal) * 1.1  ;

%yText = 1.4*valError(iConfigOptimal) ;
%yTextDelta = yText*0.36 ;
%text(xText, yText, 0, sprintf('Target Err. %2.2f',maxAllowedError))
%text(xText, yText-yTextDelta, 0, 'Optim. config.')
%text(xText, yText-1.65*yTextDelta, 0, sprintf('FD O%i /N=%i', ...
%    valOrder(iConfigOptimal),valN(iConfigOptimal)))

TITLE = sprintf('L1 Error vs Elapse Time \n Optimal FD O%i with N=%i', ...
    valOrder(iConfigOptimal),valN(iConfigOptimal))
title(TITLE, 'FontSize', 12)

ax=gca
ax.XScale='log'
ax.YScale='log'

% plot Error versus N1
subplot(2,2,1)
hold on; grid on;
xlabel('N')
ylabel('L1 Error')
title('L1 Error vs N', 'FontSize', 12)

TITLE = sprintf('L1 Error vs N \n Derivative D%i along AXIS %i', ORD, DER)
title(TITLE, 'FontSize', 12)

for ii=1:nConfig
    if strcmp(val.textdata(ii,4), 'Ac2Standard')
        if valOrder(ii) == 2
            plot(valN(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'm', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2)
        else
            plot(valN(ii), valError(ii), 'sw', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'LineWidth', 2)
        end
    else
        if valOrder(ii) == 2
            plot(valN(ii), valError(ii), '+', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN(ii), valError(ii), '+', 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN(ii), valError(ii), '+', 'MarkerEdgeColor', 'm', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN(ii), valError(ii), '+', 'MarkerEdgeColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        else
            plot(valN(ii), valError(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% plot horizontal line with allowed error
plot(maxErrX, maxErrY, '--k')

% represent optimal config with a star
plot(valN(iConfigOptimal), valError(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

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
            plot(valN(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'm', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2)
        else
            plot(valN(ii), valGpoint(ii), 'sw', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'LineWidth', 2)
        end
    else
        if valOrder(ii) == 2
            plot(valN(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'm', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        else
            plot(valN(ii), valGpoint(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% represent optimal config with a star
plot(valN(iConfigOptimal), valGpoint(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

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
            plot(valN(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'm', 'MarkerSize', 10, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2)
        else
            plot(valN(ii), valGbyte(ii), 'sw', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'LineWidth', 2)
        end
    else
        if valOrder(ii) == 2
            plot(valN(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 4
            plot(valN(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 8
            plot(valN(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'm', 'MarkerSize', 15, 'LineWidth', 2)
        elseif valOrder(ii) == 12
            plot(valN(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        else
            plot(valN(ii), valGbyte(ii), '+', 'MarkerEdgeColor', 'g', 'MarkerSize', 15, 'LineWidth', 2)
        end
    end
end

% represent optimal config with a star
plot(valN(iConfigOptimal), valGbyte(iConfigOptimal), 'pw', 'MarkerEdgeColor', 'k', 'MarkerSize', 27, 'LineWidth', 1.5)

ax=gca
ax.XScale='log'
ax.YScale='log'

% save figure
figName = sprintf('%s-Axis%i.jpg', FILE, DER) ;
print(figName, '-djpeg')


% END
