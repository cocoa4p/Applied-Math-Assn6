%INPUTS
%t: current time
%V: system state. V = [U;dUdt] where
% U and dUdt are n x 1 column vectors
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
function dVdt = string_rate_func01(t,V,string_params)
    n = string_params.n; %number of masses
    M = string_params.M; %total mass attached to the string
    Uf_func = string_params.Uf_func; %function describing motion of end point
    dUfdt_func = string_params.dUfdt_func; %time derivative of Uf
    Tf = string_params.Tf; %tension in string
    L = string_params.L; %length of string
    c = string_params.c; %damping coefficient
    dx = string_params.dx; %horizontal spacing between masses
   
    m = M/n;

    %unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);
    %compute acceleration
    U_left = [0; U(1:n-1)];
    U_right = [U(2:n); Uf];
    
    dU_left = [0; dUdt(1:n-1)];
    dU_right = [dUdt(2:n); dUfdt];

    d2Udt2 = (Tf/dx) * (U_left - 2*U + U_right)/m + (c/dx) * ((dU_left - 2*dUdt + dU_right)/m); 
    %assemble state derivative
    dVdt = [dUdt;d2Udt2];
end