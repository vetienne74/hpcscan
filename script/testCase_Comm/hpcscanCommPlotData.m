
clear all; close all ;

iSendGB_0to1 = 10 ;
iSendGPoint_0to1 = 11 ;
iSendGB_0to2 = 12 ;
iSendGPoint_0to2 = 13 ;
iSendGB_0to3 = 14 ;
iSendGPoint_0to3 = 15 ;
iSendGB_0to4 = 16 ;
iSendGPoint_0to4 = 17 ;
iSendGB_0to5 = 18 ;
iSendGPoint_0to5 = 19 ;
iSendGB_0to6 = 20 ;
iSendGPoint_0to6 = 21 ;
iSendGB_0to7 = 22 ;
iSendGPoint_0to7 = 23 ;

iSendRecvGB_0to1 = 24 ;
iSendRecvGPoint_0to1 = 25 ;
iSendRecvGB_0to2 = 26 ;
iSendRecvGPoint_0to2 = 27 ;
iSendRecvGB_0to3 = 28 ;
iSendRecvGPoint_0to3 = 29 ;
iSendRecvGB_0to4 = 30 ;
iSendRecvGPoint_0to4 = 31 ;
iSendRecvGB_0to5 = 32 ;
iSendRecvGPoint_0to5 = 33 ;
iSendRecvGB_0to6 = 34 ;
iSendRecvGPoint_0to6 = 35 ;
iSendRecvGB_0to7 = 36 ;
iSendRecvGPoint_0to7 = 37 ;

iHaloExchGB = 38 ;
iHaloExchGPoint = 39 ;

DIR  = '.' ;
% log file name with .log extension
%FILE = 'hpcscanCommShaheen' ;
FILE = 'hpcscan.perf.Comm' ;
TITLE = 'Test Case Comm' ;

pathFile = sprintf('%s/%s.log', DIR, FILE) ;
val = importdata(pathFile) ;

figure('Position',[100 100 1000 400])

% GB/s
subplot(1,2,1); hold on

yVal = [ val.data(1,iSendGB_0to1) ...
    val.data(1,iSendGB_0to2) ...
    val.data(1,iSendGB_0to3) ...
    val.data(1,iSendGB_0to4) ...
    val.data(1,iSendGB_0to5) ...
    val.data(1,iSendGB_0to6) ...
    val.data(1,iSendGB_0to7) ...
    val.data(1,iSendRecvGB_0to1) ...
    val.data(1,iSendRecvGB_0to2) ...
    val.data(1,iSendRecvGB_0to3) ...
    val.data(1,iSendRecvGB_0to4) ...
    val.data(1,iSendRecvGB_0to5) ...
    val.data(1,iSendRecvGB_0to6) ...
    val.data(1,iSendRecvGB_0to7) ...
    val.data(1,iHaloExchGB) ...
    val.data(2,iHaloExchGB) ...
    val.data(3,iHaloExchGB)
    
    ]

cVal = categorical({'Send 0->1','Send 0->2','Send 0->3','Send 0->4','Send 0->5','Send 0->6','Send 0->7'...
    'SendRecv 0<>1','SendRecv 0<>2','SendRecv 0<>3','SendRecv 0<>4','SendRecv 0<>5','SendRecv 0<>6','SendRecv 0<>7'...
    'HaloExch 1x4x2', 'HaloExch 1x2x4', 'HaloExch 2x2x2'});
cVal = reordercats(cVal,{'Send 0->1','Send 0->2','Send 0->3','Send 0->4','Send 0->5','Send 0->6','Send 0->7'...
    'SendRecv 0<>1','SendRecv 0<>2','SendRecv 0<>3','SendRecv 0<>4','SendRecv 0<>5','SendRecv 0<>6','SendRecv 0<>7'...
    'HaloExch 1x4x2', 'HaloExch 1x2x4', 'HaloExch 2x2x2'});

b = bar(cVal,yVal) ;

b.FaceColor = 'flat';
% color for SendRecv
b.CData(8,:) = [0 .5 .5] ; b.CData(9,:) = [0 .5 .5] ; b.CData(10,:) = [0 .5 .5] ; b.CData(11,:) = [0 .5 .5] ;
b.CData(12,:) = [0 .5 .5] ; b.CData(13,:) = [0 .5 .5] ; b.CData(14,:) = [0 .5 .5] ;
% color for HaloExch
b.CData(15,:) = [0.5 0.5 0.] ; b.CData(16,:) = [0.5 0.5 0.] ; b.CData(17,:) = [0.5 0.5 0.] ;

ylabel('GByte/s'); title (TITLE)
grid on

% Gpoint/s
subplot(1,2,2); hold on

yVal = [ val.data(1,iSendGPoint_0to1) ...
    val.data(1,iSendGPoint_0to2) ...
    val.data(1,iSendGPoint_0to3) ...
    val.data(1,iSendGPoint_0to4) ...
    val.data(1,iSendGPoint_0to5) ...
    val.data(1,iSendGPoint_0to6) ...
    val.data(1,iSendGPoint_0to7) ...
    val.data(1,iSendRecvGPoint_0to1) ...
    val.data(1,iSendRecvGPoint_0to2) ...
    val.data(1,iSendRecvGPoint_0to3) ...
    val.data(1,iSendRecvGPoint_0to4) ...
    val.data(1,iSendRecvGPoint_0to5) ...
    val.data(1,iSendRecvGPoint_0to6) ...
    val.data(1,iSendRecvGPoint_0to7) ...
    val.data(1,iHaloExchGPoint) ...
    val.data(2,iHaloExchGPoint) ...
    val.data(3,iHaloExchGPoint)
    
    ]

cVal = categorical({'Send 0->1','Send 0->2','Send 0->3','Send 0->4','Send 0->5','Send 0->6','Send 0->7'...
    'SendRecv 0<>1','SendRecv 0<>2','SendRecv 0<>3','SendRecv 0<>4','SendRecv 0<>5','SendRecv 0<>6','SendRecv 0<>7'...
    'HaloExch 1x4x2', 'HaloExch 1x2x4', 'HaloExch 2x2x2'});
cVal = reordercats(cVal,{'Send 0->1','Send 0->2','Send 0->3','Send 0->4','Send 0->5','Send 0->6','Send 0->7'...
    'SendRecv 0<>1','SendRecv 0<>2','SendRecv 0<>3','SendRecv 0<>4','SendRecv 0<>5','SendRecv 0<>6','SendRecv 0<>7'...
    'HaloExch 1x4x2', 'HaloExch 1x2x4', 'HaloExch 2x2x2'});

b = bar(cVal,yVal) ;

b.FaceColor = 'flat';
% color for SendRecv
b.CData(8,:) = [0 .5 .5] ; b.CData(9,:) = [0 .5 .5] ; b.CData(10,:) = [0 .5 .5] ; b.CData(11,:) = [0 .5 .5] ;
b.CData(12,:) = [0 .5 .5] ; b.CData(13,:) = [0 .5 .5] ; b.CData(14,:) = [0 .5 .5] ;
% color for HaloExch
b.CData(15,:) = [0.5 0.5 0.] ; b.CData(16,:) = [0.5 0.5 0.] ; b.CData(17,:) = [0.5 0.5 0.] ;

ylabel('GPoint/s'); title (TITLE)
grid on

figFile = sprintf('%s.jpg', FILE) ;
print('-djpeg', figFile)

% END
