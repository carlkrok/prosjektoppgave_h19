function h = plot_gaussian_2d_xz(means, C, sdwidth)
means2 = means(2);
means13 = [means(1),means(3)]';
npts=50;
axh = gca;
set(axh, 'nextplot', 'add');
% plot the gaussian fits
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig([C(1,1), C(1,3); C(3,1), C(3,3)]); 
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means13, 1, size(ap,2)); 
h = plot3(bp(1,:), ones(npts)*means2, bp(2,:), '-', 'parent', axh,'LineWidth',1.5);