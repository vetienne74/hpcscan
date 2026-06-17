
close all ;
clear all ;

DIR  = '.' ;
MARKER_SIZE = 5 ;
STAR_SIZE = 14 ;

WRITE_INFO = 0 ;
PLOT_ERROR_ENERGY = 0 ;
FIG_NAME = 'ERROR_ENERGY' ;
XMIN = 0.1 ; XMAX = 1000 ;

PLOT_ERROR_TIME   = 1 ;
FIG_NAME = 'ERROR_TIME' ;
XMIN = 1 ; XMAX = 10000 ;

maxAllowedError = 0.03 ; % rtm prod
%maxAllowedError = 0.015 ; % very accurate modelling
%maxAllowedError = 0.25 ; % rtm quick mode

figure('Position',[100 100 1100 350])

% figure composion
ifig = 1 ; arch = 1 ; % cpu cacheblk
%ifig = 2 ; arch = [1 2] ; % cpu cacheblk + gpu optim
%ifig = 3 ; arch = [1 2 3] ; % cpu cacheblk + gpu optim + nec sca
%ifig = 1 ; arch = 3 ;

% log file name with .log extension

for iarch = arch
    
    switch iarch
        
        case 1            
            FILE = 'propa.perf.icx1.2021-07-29' ;
            FILE2 = 'propa.hwCounterEmon.icx1.2021-08-04' ;
            COLORPLOT = 'b' ;
            ARCH = 'CPU' ;
            MODE = 'CacheBlk' ;
            
        case 2
            FILE = 'propa.perf.xgpv21.2021-07-26' ;
            FILE2 = 'propa.hwCounter.xgpv21.2021-07-26' ;
            COLORPLOT = '#77AC30' ;
            ARCH = 'GPU' ;
            MODE = 'Optim' ;
            
        case 3
            FILE = 'propa.perf.plconec01.2021-07-27' ;
            FILE2 = 'propa.hwCounter.plconec01.2021-07-27' ;
            COLORPLOT = 'r' ;
            ARCH = 'VE' ;
            MODE = 'NEC SCA' ;
    end
    
    pathFile = sprintf('%s/%s.log', DIR, FILE) ;
    val = importdata(pathFile) ;
    
    valTime     = val.data(:,14) ;
    valN1       = val.data(:,6) ;
    valError    = val.data(:,15) ;
    valOrder    = val.data(:,9) ;
    valGpoint   = val.data(:,12) ;
    valGbyte    = val.data(:,13) ;
    valGflop    = val.data(:,10) ;
    
    valDt       = val.data(:,17) ;
    valStableDt = val.data(:,18) ;
    
    sizeVal = size(val.data) ;
    nConfig = sizeVal(1)
    
    pathFile2 = sprintf('%s/%s.log', DIR, FILE2) ;
    val2 = importdata(pathFile2) ;
    
    if (iarch == 1)
        % log file from Emon, average perf is position 4
        % energy = aver power * elapse time
        valEnergy = val2(:,4) .* valTime / 3600 ;
    else
        valEnergy = val2.data(:,15) ;
    end
    
    %===============================================================================================
    % plot Error versus Time
    %===============================================================================================
    
    if (PLOT_ERROR_TIME)
        %subplot(2,2,[1 2])
        hold on; grid on;
        xlabel('Elapsed Time (s)')
        ylabel('L1 Error')
        title('L1 Error vs Elapsed Time', 'FontSize', 12)
        
        elapseTimeOptimal = 9.99e+9;
        
        for ii=1:nConfig
            
            % search for optimal config
            if valError(ii) < maxAllowedError
                if valTime(ii) < elapseTimeOptimal
                    iConfigOptimal = ii ;
                    elapseTimeOptimal = valTime(ii) ;
                end
            end
            
            xVal = valTime(ii) ;
            yVal = valError(ii);
            
            if strcmp(val.textdata(ii,4), 'Ac2Standard')
                if valOrder(ii) == 4
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 6
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 10
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            else
                if valOrder(ii) == 4
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 6
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 10
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            end
        end
        
        % represent optimal config with a star
        plot(valTime(iConfigOptimal), valError(iConfigOptimal), '^', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', STAR_SIZE, 'LineWidth', 1.5)
        
        % display info on best best config
        if (WRITE_INFO)
            %xText = 0.2*min(valTime(:)) ;
            xText=0.2 ;
            yText = 1.4*valError(iConfigOptimal) ;
            yTextDelta = yText*0.36 ;
            text(xText, yText, 0, sprintf('Target Err. %2.2f',maxAllowedError))
            text(xText, yText-yTextDelta, 0, 'Optim. config.')
            text(xText, yText-1.65*yTextDelta, 0, sprintf('FD O%i %s',valOrder(iConfigOptimal),...
                string(val.textdata(iConfigOptimal,4))))
            ratioCFL= valDt(iConfigOptimal) / valStableDt(iConfigOptimal) * 100;
            text(xText, yText-2.05*yTextDelta, 0, sprintf('N=%i, %2.2f%% CFL',valN1(iConfigOptimal), ratioCFL))
        end
        
        ax=gca
        ax.XScale='log'
        ax.YScale='log'
        
        axis([XMIN XMAX 1e-2 3e0])
        
        plot(xlim,[maxAllowedError maxAllowedError], '--k', 'LineWidth', 1.5)
        
    end
    
    %===============================================================================================
    % plot Energy versus Time
    %===============================================================================================
    
    if (PLOT_ERROR_ENERGY)
        %figure('Position',[100 100 1100 450])
        %subplot(2,2,[1 2])
        hold on; grid on;
        xlabel('Energy (W.h)')
        ylabel('L1 Error')
        title('L1 Error vs Energy', 'FontSize', 12)
        
        energyOptimal = 9.99e+9;
        
        totEnergyStand = 0. ;
        totEnergySplit = 0. ;
        
        for ii=1:nConfig
            
            % search for optimal config
            if valError(ii) < maxAllowedError
                if valTime(ii) < energyOptimal
                    iConfigOptimal = ii ;
                    energyOptimal = valTime(ii) ;
                end
            end
            
            xVal = valEnergy(ii) ;
            yVal = valError(ii);
            
            if strcmp(val.textdata(ii,4), 'Ac2Standard')
                if valOrder(ii) == 4
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 6
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 10
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
                totEnergyStand = totEnergyStand + xVal ;
            else
                if valOrder(ii) == 4
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 6
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 10
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
                totEnergySplit = totEnergySplit + xVal ;
            end
        end
        
        % represent optimal config with a star
        plot(valEnergy(iConfigOptimal), valError(iConfigOptimal), '^', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', STAR_SIZE, 'LineWidth', 1.5)
        
        % display info on best best config
        if (WRITE_INFO)
            %xText = 0.2*min(valTime(:)) ;
            xText=0.2 ;
            yText = 1.4*valError(iConfigOptimal) ;
            yTextDelta = yText*0.36 ;
            text(xText, yText, 0, sprintf('Target Err. %2.2f',maxAllowedError))
            text(xText, yText-yTextDelta, 0, 'Optim. config.')
            text(xText, yText-1.65*yTextDelta, 0, sprintf('FD O%i %s',valOrder(iConfigOptimal),...
                string(val.textdata(iConfigOptimal,4))))
            ratioCFL= valDt(iConfigOptimal) / valStableDt(iConfigOptimal) * 100;
            text(xText, yText-2.05*yTextDelta, 0, sprintf('N=%i, %2.2f%% CFL',valN1(iConfigOptimal), ratioCFL))
        end
        
        ax=gca
        ax.XScale='log'
        ax.YScale='log'
        
        axis([XMIN XMAX 5e-3 3e0])
        
        plot(xlim,[maxAllowedError maxAllowedError], '--k', 'LineWidth', 1.5)
                
        fprintf('%s totEnergyStand kW.h %f \n', ARCH, totEnergyStand/1000) ;
        fprintf('%s totEnergySplit kW.h %f \n', ARCH, totEnergySplit/1000) ;        
        
    end
    
    
    % plot Error versus N1
    if 0
        %subplot(2,2,1)
        hold on; grid on;
        xlabel('N')
        ylabel('L1 Error')
        title('L1 Error vs N', 'FontSize', 12)
        for ii=1:nConfig
            
            xVal = valN1(ii) ;
            yVal = valError(ii);
            
            if strcmp(val.textdata(ii,4), 'Ac2Standard')
                if valOrder(ii) == 4
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 6
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 10
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            else
                if valOrder(ii) == 4
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 6
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 10
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            end
        end
        
        %axis([500 1000 1e-3 1e0])
        
        % represent optimal config with a star
        %plot(valN1(iConfigOptimal), valError(iConfigOptimal), '^', 'MarkerFaceColor', COLORPLOT, 'MarkerEdgeColor', 'k', 'MarkerSize', STAR_SIZE, 'LineWidth', 1.5)
        
        ax=gca
        ax.XScale='log'
        ax.YScale='log'
        
    end
    
    % plot Gpoint versus N1
    if 0
        subplot(2,2,3)
        hold on; grid on;
        xlabel('N')
        ylabel('GPoint/s')
        title('GPoint/s vs N', 'FontSize', 12)
        for ii=1:nConfig
            
            xVal = valN1(ii) ;
            yVal = valGpoint(ii) ;
            
            if strcmp(val.textdata(ii,4), 'Ac2Standard')
                if valOrder(ii) == 2
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 4
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 12
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            else
                if valOrder(ii) == 2
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 4
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 12
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            end
        end
        
        % represent optimal config with a star
        plot(valN1(iConfigOptimal), valGpoint(iConfigOptimal), '^', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', STAR_SIZE, 'LineWidth', 1.5)
        
        ax=gca
        ax.XScale='log'
        %ax.YScale='log'
        
        axis([500 1000 0 32])
    end
    
    % plot Gflop versus N1
    if 0
        subplot(2,2,4)
        hold on; grid on;
        xlabel('N')
        ylabel('GFlop/s')
        title('GFlop/s vs N', 'FontSize', 12)
        for ii=1:nConfig
            
            xVal = valN1(ii) ;
            yVal = valGflop(ii) ;
            
            if strcmp(val.textdata(ii,4), 'Ac2Standard')
                if valOrder(ii) == 2
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 4
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                elseif valOrder(ii) == 12
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE+2, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            else
                if valOrder(ii) == 2
                    plot(xVal, yVal, 'd', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 4
                    plot(xVal, yVal, 's', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 8
                    plot(xVal, yVal, 'p', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                elseif valOrder(ii) == 12
                    plot(xVal, yVal, 'o', 'MarkerEdgeColor', COLORPLOT, 'MarkerFaceColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                else
                    plot(xVal, yVal, '+', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', MARKER_SIZE, 'LineWidth', 1)
                end
            end
        end
        
        % represent optimal config with a star
        plot(valN1(iConfigOptimal), valGflop(iConfigOptimal), '^', 'MarkerEdgeColor', COLORPLOT, 'MarkerSize', STAR_SIZE, 'LineWidth', 1.5)
        
        % plot mem BW
        %plot(xlim,[216 216], '--b', 'LineWidth', 2)
        %plot(xlim,[1085 1085], '--r', 'LineWidth', 2)
        
        ax=gca
        ax.XScale='log'
        %ax.YScale='log'
        
        %axis([500 1000 0 3300])
        
    end
    
end

% save figure
figName = sprintf('propaParamAnalysis%s_fig%d.jpg', FIG_NAME, ifig) ;
print(figName, '-djpeg')
