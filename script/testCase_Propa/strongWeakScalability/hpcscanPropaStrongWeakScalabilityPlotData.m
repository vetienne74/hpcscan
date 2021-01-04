
clear all; close all ;

iNproc    = 1 ;
ipropaGpointFD = 11 ;
ipropaGpointEff = 12 ; 

for ifig = 1:1
    
    if ifig == 1
        DIR  = '.' ;
        % log file name with .log extension
        %FILE = 'hpcscanPropaStrongWeakScalabilityShaheen' ;
        FILE = 'hpcscan.perf.Propa' ;
    else
        DIR  = '.' ;
        FILE = 'hpcscanPropaStrongWeakScalabilityXXX' ;
    end
    
    figure('Position',[100 100 1000 400])
    
    pathFile = sprintf('%s/%s.log', DIR, FILE) ;
    val = importdata(pathFile) ;
            
    % strong scalability
    subplot(1,2,1); hold on
    
    % O4
    xVal = val.data(1:5,iNproc) ;
    yPropaGpointEff = val.data(1:5,ipropaGpointEff) ;
    plot(xVal, yPropaGpointEff, 'ko-', 'LineWidth', 2., 'DisplayName', 'Propa O4')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal O4')
    yMax = max(yValIdeal) ;
    
    % O8
    xVal = val.data(6:10,iNproc) ;
    yPropaGpointEff = val.data(6:10,ipropaGpointEff) ;
    plot(xVal, yPropaGpointEff, 'ro-', 'LineWidth', 2., 'DisplayName', 'Propa O8')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ideal O4')
    if max(yValIdeal) > yMax 
        yMax = max(yValIdeal) ;
    end
    
    % O12
    xVal = val.data(11:15,iNproc) ;
    yPropaGpointEff = val.data(11:15,ipropaGpointEff) ;
    plot(xVal, yPropaGpointEff, 'bo-', 'LineWidth', 2., 'DisplayName', 'Propa O8')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Ideal O4')
    if max(yValIdeal) > yMax 
        yMax = max(yValIdeal) ;
    end
    
    xlabel('# MPI'); ylabel('GPoint/s'); title ('Propagator Strong Scalability')
    xMax = max(val.data(:,iNproc)) ;
    axis([1 xMax 0 yMax]) ; grid on
    legend('show','Location','northwest')
    
    % weak scalability
    subplot(1,2,2); hold on
    
    % O4
    xVal = val.data(16:20,iNproc) ;
    yPropaGpointEff = val.data(16:20,ipropaGpointEff) ;
    plot(xVal, yPropaGpointEff, 'ko-', 'LineWidth', 2., 'DisplayName', 'Propa O4')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal O4')
    
    % O8
    xVal = val.data(21:25,iNproc) ;
    yPropaGpointEff = val.data(21:25,ipropaGpointEff) ;
    plot(xVal, yPropaGpointEff, 'ro-', 'LineWidth', 2., 'DisplayName', 'Propa O8')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ideal O4')
    
    % O12
    xVal = val.data(26:30,iNproc) ;
    yPropaGpointEff = val.data(26:30,ipropaGpointEff) ;
    plot(xVal, yPropaGpointEff, 'bo-', 'LineWidth', 2., 'DisplayName', 'Propa O8')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Ideal O4')
    
    xlabel('# MPI'); ylabel('GPoint/s'); title ('Propagator Weak Scalability')
    xMax = max(val.data(:,iNproc)) ;
    axis([1 xMax 0 yMax]) ; grid on
    legend('show','Location','northwest')
    
    figFile = sprintf('%s.jpg', FILE) ;
    print('-djpeg', figFile)
    
end
