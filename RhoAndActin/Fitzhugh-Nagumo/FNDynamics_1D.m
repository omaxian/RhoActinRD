% Fitzhugh-Nagumo in 1D
rng(1); % For reproducibility
L = 20;
a = 0.5;
b = 0.5;
Iext = 0.5;
tau = 10;
Du = 0.1; % The size of the waves depends on Du
Dv = 0.01;
tf = 60;
MakeMovie = 1;

% Numerical params
dt = 0.1;
saveEvery=floor(1e-6+0.1/dt);
Nx = 400; % The grid spacing has GOT to resolve the boundary layer - scales with min(D)!
dx = L/Nx;
kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
ksq=kvals.^2;
DivFacFourier_u = (1/dt+Du*ksq);
DivFacFourier_v = (1/dt+Dv*ksq);
nSave = floor(tf/(saveEvery*dt)+1e-10);
Allu=zeros(nSave,Nx);
Allv=zeros(nSave,Nx);

% Steady states
[rts,stability,~] = findFNRoots(a,b,Iext,tau,Du,Dv,Nx);
[nRts,~]=size(rts);
rts=real(rts);
if (nRts>1)
    RtRange = max(rts)-min(rts);
    RtStart = min(rts)-1.1*RtRange;
    RtEnd = max(rts)+1.1*RtRange;
else
    RtRange = 1;
    RtStart = rts-0.5*RtRange;
    RtEnd = rts+0.5*RtRange;
end
u = RtStart(1)+(RtEnd(1)-RtStart(1)).*rand(1,Nx);
v = RtStart(2)+(RtEnd(2)-RtStart(2)).*rand(1,Nx);

if (MakeMovie)
    f=figure;
end
nSt = tf/dt;
for iT=1:nSt
    RHS_u = u-1/3*u.^3-v+Iext;
    RHS_v = 1/tau*(u+a-b*v);
    RHSHat_u = fft(u/dt+RHS_u);
    uHatNew = RHSHat_u./(DivFacFourier_u);
    u = ifft(uHatNew);
    RHSHat_v = fft(v/dt+RHS_v);
    vHatNew = RHSHat_v./(DivFacFourier_v);
    v = ifft(vHatNew);
    % Print error if it goes unstable
    if (max(abs(u(:)))>1e5)
        warning('Unstable simulation - returning')
        Stats.XCor=0;
        Stats.MeanActin=inf;
        return
    end
    if (mod(iT-1,10)==0 && MakeMovie)
        plot([(0:Nx-1) Nx]*dx,[u u(1)],[(0:Nx-1) Nx]*dx,[v v(1)]);
        hold on
        xlabel('$x$ ($\mu$m)')
        ylabel('$u,v$')
        %colorbar
        pbaspect([1 1 1])
        title(strcat('$t=$',num2str(iT*dt)))
        ylim([-3 3])
        movieframes((iT-1)/10+1)=getframe(f);
        hold off
    end
    if (mod(iT-1,saveEvery)==0)
        Allu((iT-1)/saveEvery+1,:) = u;
        Allv((iT-1)/saveEvery+1,:) = v;
    end
end

% Make kymograph
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
nexttile
imagesc(Allu)
title('$u(x,t)$')
xlabel('$x$')
ylabel('$t$')
nexttile
imagesc(Allv)
title('$v(x,t)$')
xlabel('$x$')
ylabel('$t$')