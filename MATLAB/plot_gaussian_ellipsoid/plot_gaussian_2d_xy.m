function h = plot_gaussian_2d_xy(means, C, sdwidth)
means3 = means(3);
means12 = means(1:2)';
npts=50;
axh = gca;
set(axh, 'nextplot', 'add');
% plot the gaussian fits
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(C(1:2,1:2)); 
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means12, 1, size(ap,2)); 
h = plot3(bp(1,:), bp(2,:), ones(npts)*means3, '-', 'parent', axh);