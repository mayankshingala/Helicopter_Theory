%% Calculations of Force & Moments (Integration)

ds = (R - e)/20;            % Blade is divided in 20 segments
psi = psi_all(1,:);

for i = 1 : 360

    for j = 1 : 20
        a = e + ds*(j-1);
        b = a + ds;

        for k = 1 : 6
            t(k) = ((b-a)/2)*Node(k) + ((b+a)/2);      % change of variables: new locations are t(i) from x(i)
            theta(i) = theta_0 + theta_tw*t(k) + theta_1s*sin(psi(i)) + theta_1c*cos(psi(i));   % Pilot inputs

            % Non-dimensional velocity
%             u_t(k) = t(k) + mu*sin(psi(i));
%             u_p(k) = lambda*cos(beta(i)) + t(k)*beta_star(i) + mu*cos(psi(i))*sin(beta(i));
%             u_r(k) = mu*cos(psi(i));

            % Dimensional Velocity
            UT(k) = t(k)*(omega) + mu*sin(psi(i))*(omega*R);
            UP(k) = lambda*cos(beta(i))*(omega*R) + t(k)*beta_dot(i) + mu*cos(psi(i))*sin(beta(i))*(omega*R);
            UR(k) = mu*cos(psi(i))*(omega*R);

            % dFz_dr, dFx_dr, dFr_dr
            dFz_dr(k) = 0.5*rho*C*Cla* ( ((UT(k)^2)*theta(i)) - (UP(k)*UT(k)) ) * cos(beta(i));
            dFx_dr(k) = 0.5*rho*C*Cla* ( (UP(k)*UT(k)*theta(i)) - (UP(k)^2) + (Cd*(UT(k)^2)/Cla) );
            dFr_dr(k) = -0.5*rho*C*Cla*beta(i)* ( ((UT(k)^2)*theta(i)) - (UP(k)*UT(k)) );

            % dFz, dFx, dFr (6 points for integration)
            dFz(k) = dFz_dr(k) * W(k) * ((b-a)/2);
            dFx(k) = dFx_dr(k) * W(k) * ((b-a)/2);
            dFr(k) = dFr_dr(k) * W(k) * ((b-a)/2);

            % Shear forces
            Sx_rot(k) = dFx_dr(k) * W(k) * ((b-a)/2);
            Sr_rot(k) = ( m*(omega^2)*t(k) - beta(i)*dFr_dr(k) ) * W(k) * ((b-a)/2);
            Sz_rot(k) = ( dFz_dr(k) - m*(t(k)-e)*beta_ddot(i) ) * W(k) * ((b-a)/2);

            % Blade bending loads(Rotating Frame)
%             nf_rot(k) = ((t(k)-e)*dFz_dr(k)) * W(k) * ((b-a)/2);
%             nf_rot(k) = (nu_beta^2 - 1 - (3/2)*(e/R)) * I_beta * (omega^2) * (beta(i) - beta_p);
            nf_rot(k) = ( ( (t(k)-e)*dFz_dr(k) )... 
                        - ( m*(omega^2)*t(k)*(t(k)-e)*beta(i) )....
                        - ( m* ((t(k)-e)^2) *beta_ddot(i) ) ) * W(k) * ((b-a)/2);
            nl_rot(k) = ( (t(k) - e)*dFx_dr(k) ) * W(k) * ((b-a)/2);
            nt_rot(k) = 0;      % Assume elastic axis coincide with quarter chord and pitching moment is zero

        end

        % Summation of forces and moments for one element for one azimuth (ψ) (Gaussian quadrature)
        Sx_rot_ele(j) = sum(Sx_rot);
        Sr_rot_ele(j) = sum(Sr_rot);
        Sz_rot_ele(j) = sum(Sz_rot);
        nf_rot_ele(j) = sum(nf_rot);
        nl_rot_ele(j) = sum(nl_rot);
        nt_rot_ele(j) = sum(nt_rot);

    end

    % Root shear loads in Rotating frame (for different azimuth (ψ))
    Sx(i) = sum(Sx_rot_ele);
    Sr(i) = sum(Sr_rot_ele);
    Sz(i) = sum(Sz_rot_ele);

    % Root bending loads in Rotating frame (for different azimuth (ψ))
    nf(i) = sum(nf_rot_ele);
    nl(i) = sum(nl_rot_ele);
    nt(i) = sum(nt_rot_ele);

    % Hub loads (Rotating frame)
    fx(1,i) = Sx(i);                % Values assigned to 1st blade data
    fy(1,i) = Sr(i);
    fz(1,i) = Sz(i);

    % Hub Moments (Rotating frame)
    mx(1,i) = nf(i) + e*Sz(i);
    my(1,i) = nt(i);
    mz(1,i) = -nl(i) - e*Sx(i);

end

phase_shifted_data;

% Hub loads & Moments (Fixed frame)
for blade = 1:1:4
    for l = 1 : 360             % To calculate loads over different Ψ (0° to 360°)
        psi_m(l) = psi(l) + ((blade-1)*(2*pi)/Nb);      % Azimuth for mth blade
        Fx_m_fixed(blade,l) = fy(blade,l)*cos(psi_m(l)) + fx(blade,l)*sin(psi_m(l));
        Fy_m_fixed(blade,l) = fy(blade,l)*sin(psi_m(l)) - fx(blade,l)*cos(psi_m(l));
        Fz_m_fixed(blade,l) = fz(blade,l);

        Mx_m_fixed(blade,l) = mx(blade,l)*sin(psi_m(l)) + my(blade,l)*cos(psi_m(l));
        My_m_fixed(blade,l) = -mx(blade,l)*cos(psi_m(l)) + my(blade,l)*sin(psi_m(l));
        Mz_m_fixed(blade,l) = mz(blade,l);
    end
end

% Calculation of loads 
for p = 1 : 360             % To calculate loads over different Ψ (0° to 360°)

    % Rotating Frame
    fx_plot_rot(p) = sum(fx(:,p));      % Force added for all 4 blades over different Ψ (0° to 360°)
    fy_plot_rot(p) = sum(fy(:,p));
    fz_plot_rot(p) = sum(fz(:,p));
    mx_plot_rot(p) = sum(mx(:,p));      % Moment added for all 4 blades over different Ψ (0° to 360°)
    my_plot_rot(p) = sum(my(:,p));
    mz_plot_rot(p) = sum(mz(:,p));

    % Fixed Frame
    Fx_plot_fixed(p) = sum(Fx_m_fixed(:,p));
    Fy_plot_fixed(p) = sum(Fy_m_fixed(:,p));
    Fz_plot_fixed(p) = sum(Fz_m_fixed(:,p));
    Mx_plot_fixed(p) = sum(Mx_m_fixed(:,p));
    My_plot_fixed(p) = sum(My_m_fixed(:,p));
    Mz_plot_fixed(p) = sum(Mz_m_fixed(:,p));

end