function advTerm = advectiveTerm2D(Speedx,Speedy,f,h,order)
    % Speed and f are 2D arrays

    % First do the x direction
    % Estimate speeds at midpoints
    Speedx_plus = 1/2*(Speedx+circshift(Speedx,-1,2));
    % Fluxes
    fplus_x = circshift(f,-1,2);
    if (order==1)
        Fluxes_x = (Speedx_plus > 0).*f + (Speedx_plus <=0).*fplus_x;
    elseif (order==3)
        fplus2_x = circshift(f,-2,2);
        fminus_x = circshift(f,1,2);
        PlusSide = (1/6)*(-fminus_x+5*f+2*fplus_x);
        MinusSide = (1/6)*(2*f+5*fplus_x-fplus2_x);
        Fluxes_x = (Speedx_plus > 0).*PlusSide + (Speedx_plus <=0).*MinusSide;
    end
    Fluxes_x = Fluxes_x.*Speedx_plus;
    divXTerm = (Fluxes_x - circshift(Fluxes_x,1,2))/h;

    % Now do the y direction
    % Estimate speeds at midpoints
    Speedy_plus = 1/2*(Speedy+circshift(Speedy,-1,1));
    % Fluxes
    fplus_y = circshift(f,-1,1);
    if (order==1)
        Fluxes_y = (Speedy_plus > 0).*f + (Speedy_plus <=0).*fplus_y;
    elseif (order==3)
        fplus2_y = circshift(f,-2,1);
        fminus_y = circshift(f,1,1);
        PlusSide = (1/6)*(-fminus_y+5*f+2*fplus_y);
        MinusSide = (1/6)*(2*f+5*fplus_y-fplus2_y);
        Fluxes_y = (Speedy_plus > 0).*PlusSide + (Speedy_plus <=0).*MinusSide;
    end
    Fluxes_y = Fluxes_y.*Speedy_plus;
    divYTerm = (Fluxes_y - circshift(Fluxes_y,1,1))/h;

    advTerm = divXTerm + divYTerm;
end

