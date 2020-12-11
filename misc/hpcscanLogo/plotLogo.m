
clear all ;
close all ;

FLOAT_SIZE ='float32' ;
PERC       = 1.0 ;
COLORMAP   = parula ;

NSNAP = 3 ;

figure('Position',[100 100 1000 400])

file_name  = "PropaEigenModePrn.proc0" ;
file_title = "Computed" ;

% get n1 from file .info
file_name_info = file_name + ".grid.info" ;
val = load(file_name_info) ;
N1 = val(1)
N2 = val(2)

% read data
file_name_bin = file_name + ".grid.bin" ;
f1 = fopen(file_name_bin,'r','native') ;
for isnap = 1:NSNAP
    val = fread(f1, [N1 N2], FLOAT_SIZE) ;
end
fclose(f1) ;

% plot data
%----------
hold on
%title(file_title) ;
imagesc(val)
max_val = max(max(abs(val))) ;
caxis([-PERC*max_val PERC*max_val])
axis ij; axis tight; axis off
colormap(COLORMAP); 
%colorbar ;
text(N2/2, N1/2, 0, 'hpcscan', 'FontSize', 100, 'HorizontalAlignment','center', 'Color', ...
    'k', 'FontWeight', 'bold')
colormap(COLORMAP); 

fig_name=sprintf("%s.tmp.jpg", "hpcscanLogo") ;
print(fig_name, '-djpeg') ;

% END