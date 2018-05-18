function [threshold, th_idx] = getSVDThreshold(sigma,m,n)
M         = length(sigma);
y         = sigma;
threshold = optimal_SVHT_coef(m/n,0) * median(y);
y( y < (threshold) ) = 0; th_idx = find(y==0,1,'first') ;
th_idx = th_idx-1; % last non-zero value
figure,plot(sigma,'ok'),hold on
plot([th_idx th_idx],[-20 threshold],'--r','LineWidth',2)
plot([0 th_idx],[threshold threshold],'--r','LineWidth',2)
axis([0 th_idx+100 -20 max(max([sigma]))+10])

end