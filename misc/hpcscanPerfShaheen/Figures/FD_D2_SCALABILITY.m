
clear all; close all ;

iNproc    = 1 ;
iLaplaGpointFD  = 23 ;
iLaplaGpointEff = 24 ; 

for ifig = 1:2
    
    if ifig == 1
        DIR   = '../../../script/testCase_FD_D2/results' ;
        FILE1 = 'runStrongScalaNEC.xbenchmark.perf.FD_D2.log' ;
        FILE2 = 'runWeakScalaNEC.xbenchmark.perf.FD_D2.log' ;
        FIGNAME = 'FD_D2_SCALABILITY_NEC' ;
    else
        DIR   = '../../../script/testCase_FD_D2/results' ;
        FILE1 = 'runStrongScalaShaheen.xbenchmark.perf.FD_D2.log' ;
        FILE2 = 'runWeakScalaShaheen.xbenchmark.perf.FD_D2.log' ;
        FIGNAME = 'FD_D2_SCALABILITY_SHAHEEN' ;
    end    
    
    figure('Position',[100 100 1000 400])
    
    subplot(1,2,1); hold on
    pathFile = sprintf('%s/%s', DIR, FILE1) ; 
    val = load(pathFile) ;
    xVal = val(:,iNproc) ;
    yLaplaGpointFD  = val(:,iLaplaGpointFD) ;
    yLaplaGpointEff = val(:,iLaplaGpointEff) ;
    plot(xVal, yLaplaGpointFD, 'ko-', 'LineWidth', 1., 'DisplayName','FD kernel')
    plot(xVal, yLaplaGpointEff, 'ko-', 'LineWidth', 2., 'DisplayName', 'FD + Comm.')
    yValIdeal = xVal * yLaplaGpointFD(1) ;
    plot(xVal, yValIdeal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal')
    xlabel('# MPI'); ylabel('GPoint/s'); title ('Strong scalability')
    xMax = max(xVal) ;
    yMax = max(yValIdeal) ;
    axis([1 xMax 0 yMax]) ; grid on
    legend('show','Location','northwest')
    
    subplot(1,2,2); hold on
    pathFile = sprintf('%s/%s', DIR, FILE2) ; 
    val = load(pathFile) ;
    xVal = val(:,iNproc) ;
    yLaplaGpointFD  = val(:,iLaplaGpointFD) ;
    yLaplaGpointEff = val(:,iLaplaGpointEff) ;
    plot(xVal, yLaplaGpointFD, 'ko-', 'LineWidth', 1., 'DisplayName','FD kernel')
    plot(xVal, yLaplaGpointEff, 'ko-', 'LineWidth', 2., 'DisplayName', 'FD + Comm.')
    yValIdeal = xVal * yLaplaGpointFD(1) ;   
    plot(xVal, yValIdeal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal')
    xlabel('# MPI'); ylabel('GPoint/s'); title ('Weak scalability')
    xMax = max(xVal) ;
    yMax = max(yValIdeal) ;
    axis([1 xMax 0 yMax]) ; grid on
    legend('show','Location','northwest')
    
    figFile = sprintf('%s.jpg', FIGNAME) ;
    print('-djpeg', figFile)
    
end