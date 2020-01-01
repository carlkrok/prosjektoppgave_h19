function hyz = plot_gaussian_2d_yz(means, C, sdwidth)
means1 = means(1);
means23 = means(2:3)';
npts=200;
axh = gca;
set(axh, 'nextplot', 'add');
% plot the gaussian fits
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(C(2:3,2:3)); 
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means23, 1, size(ap,2)); 
hyz = plot3(ones(npts)*means1, bp(1,:), bp(2,:), '-', 'parent', axh,'LineWidth',1.2);%,'Color','#A2142F');