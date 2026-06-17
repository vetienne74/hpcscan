
clear all; 
close all ;

PLOT_FIG = 1 ;

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

orientation1 = 'north'
orientation4 = 'south'

DIR  = '.' ;
% log file name with .log extension

figure('Position',[100 100 1200 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot n1,n2,n3,Laplacian GB/s and ratio GB/s vs Max BW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iarch = [5 6 7]
    
    switch iarch
        case 1
            FILE = 'XXX' ;
            ARCH = 'CPU' ;
            MODE = 'Baseline' ;
            COLORPLOT = 'b' ;
            WIDTHPLOT = 1.0 ;
            STYLEPLOT = '--' ;
            BW   = 281 ;
            FLOP = 9999 ;
            HLINE = '\hline' ;
            
        case 2
            FILE = 'XXX' ;
            ARCH = 'CPU' ;
            MODE = 'CacheBlk' ;
            COLORPLOT = 'b' ;
            WIDTHPLOT = 1.0 ;
            STYLEPLOT = '-' ;
            BW   = 281 ;
            FLOP = 9999 ;
            HLINE = '\hline' ;
            
        case 3
            FILE = 'XXX' ;
            ARCH = 'GPU' ;
            MODE = 'Regular' ;
            COLORPLOT = '#77AC30' ;
            WIDTHPLOT = 0.5 ;
            STYLEPLOT = '-' ;
            BW   = 900 ;
            FLOP = 9999 ;
            HLINE = '\hline' ;
            
        case 4
            FILE = 'XXX' ;
            ARCH = 'GPU' ;
            MODE = 'Optim' ;
            COLORPLOT = '#77AC30' ;
            WIDTHPLOT = 1.5 ;
            STYLEPLOT = '-' ;
            BW   = 900 ;
            FLOP = 9999 ;
            HLINE = '\hline' ;
            
        case 5
            FILE = 'hpcscanPerfFD_D2Proto27Baseline' ;
            ARCH = 'SX-Aurora'
            MODE = 'Baseline' ;           
            COLORPLOT = 'r' ;
            WIDTHPLOT = 1.0 ;
            STYLEPLOT = '--' ;
            BW   = 1351 ;
            FLOP = 4320 ;
            HLINE = '\hline' ;
            
        case 6
            FILE = 'hpcscanPerfFD_D2Proto27NEC' ;
            ARCH = 'SX-Aurora'
            MODE = 'NEC' ;
            COLORPLOT = 'r' ;
            WIDTHPLOT = 0.5 ;
            STYLEPLOT = '-' ;
            BW   = 1351 ;
            FLOP = 4320 ;
            HLINE = '\hline' ;
            
        case 7
            FILE = 'hpcscanPerfFD_D2Proto27NEC_SCA' ;
            ARCH = 'SX-Aurora'
            MODE = 'SCA' ;
            COLORPLOT = 'r' ;
            WIDTHPLOT = 1.5 ;
            STYLEPLOT = '-' ;
            BW   = 1351 ;
            FLOP = 4320 ;
            HLINE = '\hline' ;
            
    end
    
    pathFile = sprintf('%s/%s.log', DIR, FILE) ;
    val = importdata(pathFile) ;
    
    %---------------------------------------------------------------------------------
    % GB
    nOrder = 8 ;
    yD2Axis1 = val.data(1:nOrder,iD2Axis1GB) ;
    yD2Axis2 = val.data(1:nOrder,iD2Axis2GB) ;
    yD2Axis3 = val.data(1:nOrder,iD2Axis3GB) ;
    yD2Lapla = val.data(1:nOrder,iD2LaplaGB) ;
      
    i1 = 1 ; i2 = i1 + nOrder -1 ; 
    xD2Axis1 = (i1:i2) ;
    yD2(xD2Axis1) = yD2Axis1 ;
    i1 = i2 + 2 ; i2 = i1 + nOrder -1 ;
    xD2Axis2 = (i1:i2) ;
    yD2(xD2Axis2) = yD2Axis2 ;
    i1 = i2 + 2 ; i2 = i1 + nOrder -1 ;
    xD2Axis3 = (i1:i2) ;
    yD2(xD2Axis3) = yD2Axis3 ;
    i1 = i2 + 2 ; i2 = i1 + nOrder -1 ;
    xD2Lapla = (i1:i2) ;
    yD2(xD2Lapla) = yD2Lapla ;
    
    % plot mem BW
    subplot(2,1,1); hold on
  
    plot(xD2Axis1, yD2Axis1, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5)
    plot(xD2Axis2, yD2Axis2, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5)
    plot(xD2Axis3, yD2Axis3, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5)
    plot(xD2Lapla, yD2Lapla, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5)
    title ('Measured GB/s')
    grid on
    ax=gca
    ax.YScale='log'
    
    plot(xlim,[BW BW], '-.', 'Color', COLORPLOT, 'LineWidth', 0.5, 'DisplayName', ARCH)
    ylim([80 8000])
    yticks([100 200 281 500 900 1351 2000 4000 8000])
    xticks([1 2 3 4 5 6 7 8                           10 11 12 13 14 15 16 17               19 20 21 22 23 24 25 26               28 29 30 31 32 33 34 35])
    xticklabels({'O2','4','6','8','10','12','14','16','O2','4','6','8','10','12','14','16', 'O2','4','6','8','10','12','14','16', 'O2','4','6','8','10','12','14','16'})
    
    text(1,  125, 0, 'Derivative axis 1') ;
    text(10, 125, 0, 'Derivative axis 2') ;
    text(19, 125, 0, 'Derivative axis 3') ;
    text(28, 125, 0, 'Laplacian') ;
    
    % ration mem BW / peak
    subplot(2,1,2); hold on
    plot(xD2Axis1, yD2Axis1/BW, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5, 'DisplayName', strcat(ARCH, " ", MODE))
    h1=plot(xD2Axis2, yD2Axis2/BW, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5)
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2=plot(xD2Axis3, yD2Axis3/BW, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5)
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3=plot(xD2Lapla, yD2Lapla/BW, 'LineStyle', STYLEPLOT, 'Color', COLORPLOT, 'LineWidth', WIDTHPLOT, 'Marker', 'o', 'MarkerFaceColor', COLORPLOT, 'MarkerSize', 5)
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    title ('Ratio Measured/Peak GB/s')
    grid on
    ylim([0 9])
    yticks([0 1 2 3 4 5 6 7 8 9])
    xticks([1 2 3 4 5 6 7 8                           10 11 12 13 14 15 16 17               19 20 21 22 23 24 25 26               28 29 30 31 32 33 34 35])
    xticklabels({'O2','4','6','8','10','12','14','16','O2','4','6','8','10','12','14','16', 'O2','4','6','8','10','12','14','16', 'O2','4','6','8','10','12','14','16'})
    
    lgd = legend('show', 'Location', 'north', 'NumColumns', 6)
    title(lgd,'Architecture')
    lgd.FontSize = 9 ;
    
end

if PLOT_FIG
    figFile = sprintf('%s.jpg', 'testCaseFD_D2.tmp') ;
    print('-djpeg', figFile)
end

