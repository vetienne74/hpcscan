
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

for ifig = 1:1
    
    figure('Position',[100 100 1000 800])
    
    if ifig == 1
        DIR  = '.' ;
        % log file name with .log extension
        %FILE = 'hpcscanFD_D2Shaheen' ;
        FILE = 'hpcscan.perf.FD_D2' ;
        orientation1 = 'southwest' ;
        orientation2 = 'southwest' ;
        orientation3 = 'northwest' ;
        orientation4 = 'northwest' ;
    elseif ifig == 2
        DIR  = '.' ;
        FILE = 'hpcscanFD_D2Shaheen' ;
        orientation1 = 'northeast' ;
        orientation2 = 'northeast' ;
        orientation3 = 'south' ;
        orientation4 = 'south' ;
    end
    
    pathFile = sprintf('%s/%s.log', DIR, FILE) ;
    val = importdata(pathFile) ;
        
    %---------------------------------------------------------------------------------
    % GpointEff
    subplot(2,2,1); hold on    
    
    xVal = val.data(1:5,iFDOrder) ;
    
    yD2Axis1 = val.data(1:5,iD2Axis1GpointEff) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1, 'DisplayName', 'Ax1 Base')
    yD2Axis2 = val.data(1:5,iD2Axis2GpointEff) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1, 'DisplayName', 'Ax2 Base')
    yD2Axis3 = val.data(1:5,iD2Axis3GpointEff) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1, 'DisplayName', 'Ax3 Base')
    yD2Lapla = val.data(1:5,iD2LaplaGpointEff) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1, 'DisplayName', 'Lap Base')
    
    yD2Axis1 = val.data(6:10,iD2Axis1GpointEff) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 Cach')
    yD2Axis2 = val.data(6:10,iD2Axis2GpointEff) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 Cach')
    yD2Axis3 = val.data(6:10,iD2Axis3GpointEff) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 Cach')
    yD2Lapla = val.data(6:10,iD2LaplaGpointEff) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap Cach')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
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
    
    yD2Axis1 = val.data(1:5,iD2Axis1GpointFD) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1, 'DisplayName', 'Ax1 Base')
    yD2Axis2 = val.data(1:5,iD2Axis2GpointFD) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1, 'DisplayName', 'Ax2 Base')
    yD2Axis3 = val.data(1:5,iD2Axis3GpointFD) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1, 'DisplayName', 'Ax3 Base')
    yD2Lapla = val.data(1:5,iD2LaplaGpointFD) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1, 'DisplayName', 'Lap Base')
    
    yD2Axis1 = val.data(6:10,iD2Axis1GpointFD) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 Cach')
    yD2Axis2 = val.data(6:10,iD2Axis2GpointFD) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 Cach')
    yD2Axis3 = val.data(6:10,iD2Axis3GpointFD) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 Cach')
    yD2Lapla = val.data(6:10,iD2LaplaGpointFD) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap Cach')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
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
    
    yD2Axis1 = val.data(1:5,iD2Axis1Gflop) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1, 'DisplayName', 'Ax1 Base')
    yD2Axis2 = val.data(1:5,iD2Axis2Gflop) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1, 'DisplayName', 'Ax2 Base')
    yD2Axis3 = val.data(1:5,iD2Axis3Gflop) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1, 'DisplayName', 'Ax3 Base')
    yD2Lapla = val.data(1:5,iD2LaplaGflop) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1, 'DisplayName', 'Lap Base')
    
    yD2Axis1 = val.data(6:10,iD2Axis1Gflop) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 Cach')
    yD2Axis2 = val.data(6:10,iD2Axis2Gflop) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 Cach')
    yD2Axis3 = val.data(6:10,iD2Axis3Gflop) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 Cach')
    yD2Lapla = val.data(6:10,iD2LaplaGflop) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap Cach')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
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
    
    yD2Axis1 = val.data(1:5,iD2Axis1GB) ;
    plot(xVal, yD2Axis1, 'ko-', 'LineWidth', 1, 'DisplayName', 'Ax1 Base')
    yD2Axis2 = val.data(1:5,iD2Axis2GB) ;
    plot(xVal, yD2Axis2, 'bo-', 'LineWidth', 1, 'DisplayName', 'Ax2 Base')
    yD2Axis3 = val.data(1:5,iD2Axis3GB) ;
    plot(xVal, yD2Axis3, 'go-', 'LineWidth', 1, 'DisplayName', 'Ax3 Base')
    yD2Lapla = val.data(1:5,iD2LaplaGB) ;
    plot(xVal, yD2Lapla, 'mo-', 'LineWidth', 1, 'DisplayName', 'Lap Base')
    
    yD2Axis1 = val.data(6:10,iD2Axis1GB) ;
    plot(xVal, yD2Axis1, 'ko--', 'LineWidth', 1.5, 'DisplayName', 'Ax1 Cach')
    yD2Axis2 = val.data(6:10,iD2Axis2GB) ;
    plot(xVal, yD2Axis2, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Ax2 Cach')
    yD2Axis3 = val.data(6:10,iD2Axis3GB) ;
    plot(xVal, yD2Axis3, 'go--', 'LineWidth', 1.5, 'DisplayName', 'Ax3 Cach')
    yD2Lapla = val.data(6:10,iD2LaplaGB) ;
    plot(xVal, yD2Lapla, 'mo--', 'LineWidth', 1.5, 'DisplayName', 'Lap Cach')
    
    minY = 0 ; maxY = 0 ; minX = 0 ;
    maxX = max(xVal) ;
    maxY = max(maxY, max(yD2Axis1)) ;
    maxY = max(maxY, max(yD2Axis2)) ;
    maxY = max(maxY, max(yD2Axis3)) ;
    maxY = max(maxY, max(yD2Lapla)) ;    
    
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
    figFile = sprintf('%s.jpg', FILE) ;
    print('-djpeg', figFile)
    
end
