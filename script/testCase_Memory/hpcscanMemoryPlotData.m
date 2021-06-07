
clear all; close all ;

iNthread = 5 ;
iFillGridGB = 10 ;
iFillGridGpoint = 11 ;
iCopyGridGB = 12 ;
iCopyGridGpoint = 13 ;
iAddGridGB = 14 ;
iAddGridGpoint = 15 ;
iMultiplyGridGB = 16 ;
iMultiplyGridGpoint = 17 ;
iAddUpdateGridGB = 18 ;
iAddUpdateGridGpoint = 19 ;

DIR  = '.' ;
% log file name with .log extension
%FILE = 'hpcscanMemoryShaheen' ;
FILE = 'hpcscan.perf.Memory' ;
TITLE = 'Test Case Memory / scalability ' ;

pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

figure('Position',[100 100 1000 400])

subplot(1,2,1); hold on

xVal = val.data(:,iNthread) ;
plot(xVal, val.data(:,iFillGridGB), 'ko-', 'LineWidth', 1., 'DisplayName','Fill')
plot(xVal, val.data(:,iCopyGridGB), 'bo-', 'LineWidth', 1., 'DisplayName','Copy')
plot(xVal, val.data(:,iAddGridGB), 'ro-', 'LineWidth', 1., 'DisplayName','Add')
plot(xVal, val.data(:,iMultiplyGridGB) , 'go--', 'LineWidth', 1., 'DisplayName','Mul')
plot(xVal, val.data(:,iAddUpdateGridGB) , 'mo-', 'LineWidth', 1., 'DisplayName','AddUpd')
xlabel('# threads'); ylabel('GByte/s'); title (TITLE)
grid on
legend('show','Location','northwest')

subplot(1,2,2); hold on

xVal = val.data(:,iNthread) ;
plot(xVal, val.data(:,iFillGridGpoint), 'ko-', 'LineWidth', 1., 'DisplayName','Fill')
plot(xVal, val.data(:,iCopyGridGpoint), 'bo-', 'LineWidth', 1., 'DisplayName','Copy')
plot(xVal, val.data(:,iAddGridGpoint) , 'ro-', 'LineWidth', 1., 'DisplayName','Add')
plot(xVal, val.data(:,iMultiplyGridGpoint) , 'go--', 'LineWidth', 1., 'DisplayName','Mul')
plot(xVal, val.data(:,iAddUpdateGridGpoint) , 'mo-', 'LineWidth', 1., 'DisplayName','AddUpd')
xlabel('# threads'); ylabel('GPoint/s'); title (TITLE)
grid on
legend('show','Location','northwest')

figFile = sprintf('%s.jpg', FILE) ;
print('-djpeg', figFile)

% END
