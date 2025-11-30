function string_simulation_template01()
    num_masses = 2;
    total_mass = 5;
    tension_force = 5;
    string_length = 10;
    damping_coeff = .001;

    dx = string_length/(num_masses+1);
    
    amplitude_Uf = 3;
    omega_Uf = 6;



    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    
    %initial conditions
    U0 = zeros(num_masses,1); % initial disp
    dUdt0 = zeros(num_masses,1); % velocity
    V0 = [U0;dUdt0];
    tspan = linspace(0,20,2000);

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
    figure;
    h = plot(xlist, [0; Ulist(1,:)'; Uf_func(tlist(1))], 'LineWidth', 2);
    ylim([-5 5]);
    xlabel('x'); ylabel('displacement');
    
    for k = 1:length(tlist)
        uL = 0;
        uInterior = Ulist(k,:)';
        uR = Uf_func(tlist(k));
    
        y_now = [uL; uInterior; uR];
    
        set(h, 'YData', y_now);
        drawnow;
    end
end