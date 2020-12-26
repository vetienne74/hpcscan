
clear all; close all ;

iNproc    = 1 ;
ipropaGpointFD = 11 ;
ipropaGpointEff = 12 ; 

for ifig = 1:1
    
    if ifig == 1
        DIR  = '.' ;
        FILE = 'hpcscanPropaStrongWeakScalabilityShaheen' ;
    else
        DIR  = '.' ;
        FILE = 'hpcscanPropaStrongWeakScalabilityXXX' ;
    end
    
    figure('Position',[100 100 1000 400])
    
    pathFile = sprintf('%s/%s.log', DIR, FILE) ;
    val = importdata(pathFile) ;
            
    % strong scalability
    subplot(1,2,1); hold on
    xVal = val.data(1:5,iNproc) ;
    yPropaGpointFD  = val.data(1:5,ipropaGpointFD) ;
    yPropaGpointEff = val.data(1:5,ipropaGpointEff) ;
    %plot(xVal, yPropaGpointFD, 'ko-', 'LineWidth', 1., 'DisplayName','FD kernel')
    plot(xVal, yPropaGpointEff, 'ko-', 'LineWidth', 2., 'DisplayName', 'Propagator')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal')
    xlabel('# MPI'); ylabel('GPoint/s'); title ('Strong scalability')
    xMax = max(xVal) ;
    yMax = max(yValIdeal) ;
    axis([1 xMax 0 yMax]) ; grid on
    legend('show','Location','northwest')
    
    % weak scalability
    subplot(1,2,2); hold on
    xVal = val.data(6:10,iNproc) ;
    yPropaGpointFD  = val.data(6:10,ipropaGpointFD) ;
    yPropaGpointEff = val.data(6:10,ipropaGpointEff) ;
    %plot(xVal, yPropaGpointFD, 'ko-', 'LineWidth', 1., 'DisplayName','FD kernel')
    plot(xVal, yPropaGpointEff, 'ko-', 'LineWidth', 2., 'DisplayName', 'Propagator')
    yValIdeal = xVal * yPropaGpointEff(1) ;
    plot(xVal, yValIdeal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal')
    xlabel('# MPI'); ylabel('GPoint/s'); title ('Weak scalability')
    xMax = max(xVal) ;
    yMax = max(yValIdeal) ;
    axis([1 xMax 0 yMax]) ; grid on
    legend('show','Location','northwest')
    
    figFile = sprintf('%s.jpg', FILE) ;
    print('-djpeg', figFile)
    
end