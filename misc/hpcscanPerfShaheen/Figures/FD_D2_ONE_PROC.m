
clear all; close all ;

iFDOrder          = 9 ;
iD2Axis1Gflop     = 10 ;
iD2Axis1GpointFD  = 11 ;
iD2Axis1GpointEff = 12 ;
iD2Axis1GB        = 13 ;

iD2Axis2Gflop     = 14 ;
iD2Axis2GpointFD  = 15 ;
iD2Axis2GpointEff = 16 ;
iD2Axis2GB        = 17 ;

iD2Axis3Gflop     = 18 ;
iD2Axis3GpointFD  = 19 ;
iD2Axis3GpointEff = 20 ;
iD2Axis3GB        = 21 ;

iD2LaplaGflop     = 22 ;
iD2LaplaGpointFD  = 23 ;
iD2LaplaGpointEff = 24 ;
iD2LaplaGB        = 25 ;

for ifig = 1:2
    
    figure('Position',[100 100 1000 800])
    
    if ifig == 1
        DIR   = '../../../script/testCase_FD_D2/results' ;
        %FILE1 = 'runSmallGridNEC.xbenchmark.perf.FD_D2.log' ;        
        %FILE2 = 'runMediumGridNEC.xbenchmark.perf.FD_D2.log' ;
        FILE1 = 'runSmallGridNEC.Baseline.xbenchmark.perf.FD_D2.log' ;        
        FILE2 = 'runMediumGridNEC.Baseline.xbenchmark.perf.FD_D2.log' ;
        orientation1 = 'northeast' ;
        orientation2 = 'northeast' ;
        orientation3 = 'south' ;
        orientation4 = 'south' ;
        FIGNAME = 'FD_D2_ONE_PROC_NEC' ;
    elseif ifig == 2
        DIR   = '../../../script/testCase_FD_D2/results' ;
        FILE1 = 'runSmallGridNEC.NEC_SCA.xbenchmark.perf.FD_D2.log' ;        
        FILE2 = 'runMediumGridNEC.NEC_SCA.xbenchmark.perf.FD_D2.log' ;
        orientation1 = 'northeast' ;
        orientation2 = 'northeast' ;
        orientation3 = 'south' ;
        orientation4 = 'south' ;
        FIGNAME = 'FD_D2_ONE_PROC_NEC_SCA' ;
    elseif ifig == 3
        DIR   = '../../../script/testCase_FD_D2/results' ;
        FILE1 = 'runSmallGridShaheen.xbenchmark.perf.FD_D2.log' ;
        FILE2 = 'runMediumGridShaheen.xbenchmark.perf.FD_D2.log' ;
        orientation1 = 'east' ;
        orientation2 = 'east' ;    
        orientation3 = 'northwest' ;
        orientation4 = 'northwest' ;
        FIGNAME = 'FD_D2_ONE_PROC_SHAHEEN' ;
    elseif ifig == 4
        DIR   = '../../../script/testCase_FD_D2/results' ;
        FILE1 = 'runSmallGridShaheen.xbenchmark.perf.FD_D2_CacheBlk.log' ;
        FILE2 = 'runMediumGridShaheen.xbenchmark.perf.FD_D2_CacheBlk.log' ;
        orientation1 = 'southwest' ;
        orientation2 = 'southwest' ;    
        orientation3 = 'northwest' ;
        orientation4 = 'northwest' ;       
        FIGNAME = 'FD_D2_ONE_PROC_SHAHEEN_CACHEBLK' ;
    end
    
    %---------------------------------------------------------------------------------
    % GpointEff
    subplot(2,2,1); hold on    
    
    pathFile = sprintf('%s/%s', DIR, FILE1) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1GpointEff) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Ax1 S')
    yD2Axis2 = val(:,iD2Axis2GpointEff) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Ax2 S')
    yD2Axis3 = val(:,iD2Axis3GpointEff) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1.5, 'DisplayName', 'Ax3 S')
    yD2Lapla = val(:,iD2LaplaGpointEff) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1.5, 'DisplayName', 'Lap S')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
    pathFile = sprintf('%s/%s', DIR, FILE2) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1GpointEff) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 M')
    yD2Axis2 = val(:,iD2Axis2GpointEff) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 M')
    yD2Axis3 = val(:,iD2Axis3GpointEff) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 M')
    yD2Lapla = val(:,iD2LaplaGpointEff) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap M')
    
    xlabel('FD order'); ylabel('GPoint/s'); title ('GPoint/s Eff.')
    grid on    
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ; 
    axis([minX, maxX*1.1, minY, maxY*1.1])
    
    lgd = legend('show', 'Location', orientation1)
    lgd.FontSize = 7 ;
    
    %---------------------------------------------------------------------------------
    % GpointFD
    subplot(2,2,2); hold on    
    
    pathFile = sprintf('%s/%s', DIR, FILE1) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1GpointFD) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Ax1 S')
    yD2Axis2 = val(:,iD2Axis2GpointFD) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Ax2 S')
    yD2Axis3 = val(:,iD2Axis3GpointFD) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1.5, 'DisplayName', 'Ax3 S')
    yD2Lapla = val(:,iD2LaplaGpointFD) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1.5, 'DisplayName', 'Lap S')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
    pathFile = sprintf('%s/%s', DIR, FILE2) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1GpointFD) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 M')
    yD2Axis2 = val(:,iD2Axis2GpointFD) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 M')
    yD2Axis3 = val(:,iD2Axis3GpointFD) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 M')
    yD2Lapla = val(:,iD2LaplaGpointFD) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap M')
    
    xlabel('FD order'); ylabel('GPoint/s'); title ('GPoint/s FD')
    grid on    
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ; 
    axis([minX, maxX*1.1, minY, maxY*1.1])
    
    lgd = legend('show', 'Location', orientation2)
    lgd.FontSize = 7 ;

    %---------------------------------------------------------------------------------
    % GfLop
    subplot(2,2,3); hold on    
    
    pathFile = sprintf('%s/%s', DIR, FILE1) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1Gflop) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Ax1 S')
    yD2Axis2 = val(:,iD2Axis2Gflop) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Ax2 S')
    yD2Axis3 = val(:,iD2Axis3Gflop) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1.5, 'DisplayName', 'Ax3 S')
    yD2Lapla = val(:,iD2LaplaGflop) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1.5, 'DisplayName', 'Lap S')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
    pathFile = sprintf('%s/%s', DIR, FILE2) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1Gflop) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 M')
    yD2Axis2 = val(:,iD2Axis2Gflop) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 M')
    yD2Axis3 = val(:,iD2Axis3Gflop) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 M')
    yD2Lapla = val(:,iD2LaplaGflop) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap M')
    
    xlabel('FD order'); ylabel('GFLOP/s'); title ('GFLOP/s')
    grid on    
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ; 
    axis([minX, maxX*1.1, minY, maxY*1.1])
    
    lgd = legend('show', 'Location', orientation3)
    lgd.FontSize = 7 ;
    
    %---------------------------------------------------------------------------------
    % GB
    subplot(2,2,4); hold on    
    
    pathFile = sprintf('%s/%s', DIR, FILE1) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1GB) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Ax1 S')
    yD2Axis2 = val(:,iD2Axis2GB) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Ax2 S')
    yD2Axis3 = val(:,iD2Axis3GB) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1.5, 'DisplayName', 'Ax3 S')
    yD2Lapla = val(:,iD2LaplaGB) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1.5, 'DisplayName', 'Lap S')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
    pathFile = sprintf('%s/%s', DIR, FILE2) ;
    val = load(pathFile) ;
    xVal = val(:,iFDOrder) ;
    yD2Axis1 = val(:,iD2Axis1GB) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 M')
    yD2Axis2 = val(:,iD2Axis2GB) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 M')
    yD2Axis3 = val(:,iD2Axis3GB) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 M')
    yD2Lapla = val(:,iD2LaplaGB) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap M')
    
    xlabel('FD order'); ylabel('GB/s'); title ('GB/s')
    grid on    
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ; 
    axis([minX, maxX*1.1, minY, maxY*1.1])
    
    lgd = legend('show', 'Location', orientation4)
    lgd.FontSize = 7 ;
    
    %---------------------------------------------------------------------------------
    figFile = sprintf('%s.jpg', FIGNAME) ;
    print('-djpeg', figFile)
    
end