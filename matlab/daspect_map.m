function [r] = daspect_map

lats=get(gca,'YLim');
lons=get(gca,'XLim');
r=cos(pi*mean(lats)/180);

daspect([1 r 1]);
pbaspect([1 r 1]);
xlim(lons);
ylim(lats);
