% 2D heat equation using Fourier method with exact integration
% in Fourier space
L=10;
D=0.1;
gw = 0.1;
dx=gw/2;
Nx=L/dx; % The grid spacing 
x=(0:Nx-1)*dx;
y=(0:Nx-1)*dx;
[xg,yg]=meshgrid(x,y);
kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
[kx,ky]=meshgrid(kvals);
ksq=kx.^2+ky.^2;

tf=2;
u0 = sin(4*pi*xg/L).*cos(6*pi*yg/L);
SqDToSrc = (xg-L/2).^2+(yg-L/2).^2;
f0=1;
f = f0/(2*pi*gw^2)*exp(-SqDToSrc/(2*gw^2));

u=u0;
u0hat = fft2(u0);
% Fourier method (first order in time)
dt=0.001;
nSt = tf/dt;
for iT=1:nSt
    t=(iT-1)*dt;
    fhat = fft2(f);
    uhat = fft2(u);
    uhat = exp(-D*ksq*dt).*(uhat + dt*fhat);
    u = ifft2(uhat);
end
uichat = u0hat.*exp(-D*ksq*tf);
uicReal = ifft2(uichat);
uhat2 = fhat.*(-exp(-D*ksq*tf) + 1)./(D*ksq);
uhat2(1,1)=fhat(1,1)*tf;
u2 = ifft2(uhat2)+uicReal;
max(abs(u(:)-u2(:)))
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
imagesc(u0)
colorbar
nexttile
imagesc(u)
colorbar

% Fundamental solution (here is where you can coarsen the grid)
Gg0 = f0*exp(-SqDToSrc/(2*gw^2+4*tf*D))/(4*pi*tf*D+2*pi*gw^2);
uFund = -f0/(4*D*pi)*(expint(SqDToSrc/(2*gw^2))-expint(SqDToSrc/(2*gw^2+4*D*tf)));
uFund(SqDToSrc==0)=log(1 + (2*D*tf)/gw^2)/(4*D*pi);
uFund = uFund + uicReal;
nexttile
imagesc(uFund)
colorbar
max(abs(uFund(:)-u2(:)))