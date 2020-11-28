
close all ;
clear all ;

DATA = 'runMars'
%DATA = 'runNEC'
FILE = sprintf('%s.out', DATA) ;
val = importdata(FILE) ;

valTime   = val.data(:,14) ;
valN1     = val.data(:,6) ;
valError  = val.data(:,15) ;
valOrder  = val.data(:,9) ;
valGpoint = val.data(:,12) ;
valGbyte  = val.data(:,13) ;

sizeVal = size(val.data) ;
nConfig = sizeVal(1)

figure('Position',[100 100 1100 900])

% plot Error versus N1
subplot(2,2,1)
hold on; grid on;
xlabel('N1')
ylabel('Error')
title('Error vs N1', 'FontSize', 12)
for ii=1:nConfig
    if strcmp(val.textdata(ii,2), 'CacheBlk')
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
ax=gca
ax.XScale='log'
ax.YScale='log'

% plot Error versus Time
subplot(2,2,2)
hold on; grid on;
xlabel('Time (s)')
ylabel('Error')
title('Error vs Compute Time', 'FontSize', 12)
for ii=1:nConfig
    if strcmp(val.textdata(ii,2), 'CacheBlk')
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
ax=gca
ax.XScale='log'
ax.YScale='log'

% plot Gpoint versus N1
subplot(2,2,3)
hold on; grid on;
xlabel('N1')
ylabel('GPoint/s')
title('GPoint/s vs N1', 'FontSize', 12)
for ii=1:nConfig
    if strcmp(val.textdata(ii,2), 'CacheBlk')
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
ax=gca
ax.XScale='log'
ax.YScale='log'

% plot Gbyte versus N1
subplot(2,2,4)
hold on; grid on;
xlabel('N1')
ylabel('GByte/s')
title('GByte/s vs N1', 'FontSize', 12)
for ii=1:nConfig
    if strcmp(val.textdata(ii,2), 'CacheBlk')
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
ax=gca
ax.XScale='log'
ax.YScale='log'

% save figure
figName = sprintf('%s_AccPerf.tmp.jpg', DATA) ;
print(figName, '-djpeg')


% END