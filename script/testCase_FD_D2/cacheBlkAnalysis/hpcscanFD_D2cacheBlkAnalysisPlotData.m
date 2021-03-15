
clear all; close all ;

%--------------------------------------------------------------------------
% start input parameters
DIR  = '.' ;

% log file name with .log extension
%FILE = 'hpcscanFD_D2Shaheen' ;
FILE = 'hpcscan.perf.FD_D2' ;

orientation1 = 'southwest' ;
orientation2 = 'southwest' ;
orientation3 = 'northwest' ;
orientation4 = 'northwest' ;

% indicate below the range of parameters in the log file
fdOrderInLog = [4,8] ;
cb1InLog     = 9999 ;
cb2InLog     = 1:10 ;
cb3InLog     = 1:10 ;
% end input parameters
%--------------------------------------------------------------------------

% allocate tables to store perf.
nfdOrderInLog = numel(fdOrderInLog)
ncb1InLog     = numel(cb1InLog)
ncb2InLog     = numel(cb2InLog)
ncb3InLog     = numel(cb3InLog)
perfMapGpoint = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;
perfMapGflop  = zeros(nfdOrderInLog, ncb1InLog, ncb2InLog, ncb3InLog) ;

% index of parameters in log file
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

iCB1              = 26 ;
iCB2              = 27 ;
iCB3              = 28 ;

% load table
pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

% build perf maps
sizeLog = size(val.data)
if (sizeLog(1) ~= numel(perfMapGpoint))
    disp('*** ERROR, LOG FILE INCONSISTENT WITH PARAM DESCRIPTION ! ***')
    return
else
    disp('START BUILD PERF MAPS...')
end

for iconfig = 1:sizeLog(1)
    % get index in perf map
    idxFdOrder = find(fdOrderInLog(:) == val.data(iconfig,iFDOrder)) 
    idxCb1     = find(cb1InLog(:) == val.data(iconfig,iCB1)) 
    idxCb2     = find(cb2InLog(:) == val.data(iconfig,iCB2)) 
    idxCb3     = find(cb3InLog(:) == val.data(iconfig,iCB3)) 
    
    % store GFlop
    perfMapGflop(idxFdOrder, idxCb1, idxCb2, idxCb3) = val.data(iconfig,iD2LaplaGflop) ;
    
    % store GPoint
end

figure
perfMap = reshape(perfMapGflop(1, 1, :, :), [ncb2InLog ncb3InLog])
surf(perfMap)
view(0,90)
axis tight
colormap(jet)

return

fdOrderVal = val.data(1:5,iFDOrder) ;
cb1Val     = val.data(1:5,iCB1) ;
cb2Val     = val.data(1:5,iCB2) ;
cb3Val     = val.data(1:5,iCB3) ;

return

figure('Position',[100 100 1000 800])

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

%end
