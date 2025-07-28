% Improved continuum model with advection in tau direction
L=20;
Nx=40; 
advorder=3;

dt = 1e-1; % Stability limit is 1
tf = 10;
saveEvery=floor(1e-6+1/dt);

% Initialize grid
dx=L/Nx;
x1=(0:Nx-1)*dx;
y1=(0:Nx-1)*dx;
[x,y]=meshgrid(x1,y1);
Speedx = sin(6*pi*y/L);
Speedy = cos(4*pi*x/L);
u = cos(2*pi*x/L).*sin(4*pi*y/L);
kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
[kx,ky]=meshgrid(kvals);
ksq=kx.^2+ky.^2;
% Zero unpaired mode
kx(:,Nx/2+1)=0;
ky(Nx/2+1,:)=0;
nSt = floor(tf/dt+1e-10);

% Run simulation
for iT=1:nSt
    t = (iT-1)*dt;
    RHS = ((L + 4*pi*cos((4*pi*x)/L)).*cos(t - (2*pi*x)/L).*...
        cos(t + (4*pi*y)/L) + sin(t - (2*pi*x)/L).* ...
        (-L + 2*pi*sin((6*pi*y)/L)).*sin(t + (4*pi*y)/L))/L;
    advTerm = advectiveTerm2D(Speedx,Speedy,u,dx,advorder);
    u = u + dt*(-advTerm+RHS);
    % imagesc(x(1,:),y(:,1)',u)
    % drawnow
end
% Compare to exact solution
uex = cos(2*pi*x/L-t).*sin(4*pi*y/L+t);
%uex=cos(2*pi*x/L).*sin(4*pi*y/L);
max(abs(u(:)-uex(:)))