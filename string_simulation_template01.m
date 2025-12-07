function string_simulation_template01()
    num_masses = 8;
    total_mass = 8;
    tension_force = 5;
    string_length = 10;
    damping_coeff = .0001;

    dx = string_length/(num_masses+1);
    
    

    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);

    [V,D] = eig(K_mat,M_mat);

    omega_index = 5;

    mode_shape = V(:,omega_index);

    omega_n = sqrt(-D(omega_index,omega_index));
    
    
    amplitude_Uf = 3;
    omega_Uf = omega_n;

    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);

    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);

    %initial conditions
    U0 = zeros(num_masses,1); % initial disp
    dUdt0 = zeros(num_masses,1); % velocity
    V0 = [U0;dUdt0];
    tspan = [0,30*2*pi/omega_Uf];

    %run the integration
    [tlist,Vlist] = ode45(my_rate_func,tspan,V0); 
    %  input all  time values i want for ode 45
    %your code to generate an animation of the system
    Ulist = Vlist(:, 1:num_masses); % displacements
    %initialize plot
    % call plot and  save as output of  variale
    %  for loop that goes thorugh everytime step and updates plot for each
    %  time  setpggtf
    %  set  x data and  y data

    ymax = max(max(Ulist));
    ymin = min(min(Ulist));

    yq = max([abs(ymax),abs(ymin),amplitude_Uf]);

    figure;
    hold on

    mode_factor = max(abs(ymax),abs(ymin)) / max(abs(mode_shape));
    plot(xlist, mode_factor*[0; mode_shape; 0],'o-', 'LineWidth', 2, ...
        'markerfacecolor','g','markeredgecolor','g','Color','b','MarkerSize',5);

    h = plot(xlist, [0; Ulist(1,:)'; Uf_func(tlist(1))],'o-', 'LineWidth', 2, ...
        'markerfacecolor','r','markeredgecolor','r','Color','k','MarkerSize',5);
    title('Predicted Mode Shape vs Actual Vibration');
    legend(sprintf('Masses %d: Mode = %d', num_masses, omega_index));

    

    ylim(1.2*[-yq, yq]);
    xlim([0,max(xlist)]);
    xlabel('x'); ylabel('displacement');
    
    
    dT = .3/omega_Uf;
    t_index = 1;
    k = 1;

    while t_index <=  length(tlist)
        
        uL = 0;
        uInterior = Ulist(t_index,:)';
        uR = Uf_func(tlist(t_index));
    
        y_now = [uL; uInterior; uR];
    
        set(h, 'YData', y_now);
        drawnow;

        k = k+1;

        while t_index <=  length(tlist) &&  tlist(t_index)<k*dT
            t_index = t_index + 1;
        end

        pause(.03);
    end
end