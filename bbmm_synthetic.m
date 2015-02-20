% Test the Brownian Bridges Movement Model with synthetic data 

d.x = [0,0; 0.5,0.6; 1,1; 0.2,0.5; 0,0; 0.5,0.1; 1,0]; % coordinates
d.t = [0;   1;       2;   3;       4;   5;       6]; % observation times
d.s = 0.001;        % observation error (can be global or per-sample)

p_xt = brownianb(d);

x = linspace(min(d.x(:,1))-0.2, max(d.x(:,1))+0.2, 100);
y = linspace(min(d.x(:,2))-0.2, max(d.x(:,2))+0.2, 100);
[xx,yy] = meshgrid(x,y);

% expected occupation time is the integral over time of
% the probability of the particle being in a given point

tlim = [min(d.t), max(d.t)]; 
p = arrayfun(@(x,y) integral(@(t)p_xt([x,y],t), tlim(1), tlim(2)), xx, yy);
p = (p-min(p(:))) / (max(p(:))-min(p(:)));

%%
disp('plot')

figure;
imagesc(y,x,log(p)');
axis xy
drawnow

%%

p = (p-min(p(:))) / (max(p(:))-min(p(:)));
m = [xx(:), yy(:), p(:)];
m = m(p(:) > 0.3, :);

% Output file, draggable to the web interface
dlmwrite('web/data/demo0.csv', m, 'precision', '%.8f');
