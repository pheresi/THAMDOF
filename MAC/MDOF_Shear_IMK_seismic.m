function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, ug_up] = MDOF_Shear_IMK_seismic(h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck)
% [u, v, a_t, fs_st, fs, T, phi, M, K, C, time] = MDOF_ShearBiLinear_seismic(h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck)
% 
% Computes the seismic response of a nonlinear Ibarra-Medina-Krawinkler 
% MDOF system, using Newmark's method, with Newton-Raphson iterations. 
% (Chopra 4th Edition, Table 16.3.3). The MDOF system is considered as a 
% shear building. It uses constant (average) acceleration method.
% 
% Inputs:
%   - h : Vector with story heights, in [length].
%   - wi: Vector with story weights, from 1 to the roof (end of the
%         vector), in [force].
%   - Pi: Vector with P forces for considering P-delta effects, in [force].
%         If Pi = 0, no P-delta effects are considered.
%   - k: Vector with story stiffnesses, from 1 to the roof (end of the
%        vector), in [force/length].
%   - Xi: Vector with damping ratio of each mode.
%   - ug: Row vector containing the acceleration time history, in [g].
%   - dt: Time step. Recommendation: Use dt = min[dt_crit,T/10,sampling of ug]
%           Where   dt_crit = inf ; for Average Acceleration method
%                   dt_crir = 0.551*min(T) ; for Linear Acceleration method
%   - do: Vector with initial displacement, in [length].
%   - Vo: Vector with initial velocity, in [length/sec].
%   - Fy: Vector of Yield Force of each story i, in [force].
%   - a_s: Vector of hardening stiffness ratio of each story.
%   - dcdy: Vector of ductility capacity (dc/dy) for capping point of each story.
%   - a_c: Vector of post-capping stiffness ratio of each story.
%   - tol: Tolerance of the Newton-Raphson iterations. 
%           Recommendation: Use tol = 1E-3 to 1E-8.
%   - g: Standard gravity acceleration, in [length/sec^2].
%   - Name: Name of the ground motion [string].
%   - StrengthLimitCheck:   1 if strength limit is considered (recommended)
%                           0 if strength limit is NOT considered
% 
% Outputs:
%   - u: Matrix with Relative Displacement time history, in [length].
%   - v: Matrix with Relative Velocity time history, in [length/sec].
%   - a_t: Matrix with Absolute Acceleration time history, in [g].
%   - fs_st: Matrix with Story Restoring Force (Shear Force) time history, in [force].
%   - fs: Matrix with Story Forces time history, in [force].
%   - T: Vector of periods of the structure, in [sec].
%   - phi: Matrix with modal shapes of the structure.
%   - M, K, C: Mass, Stiffness and Damping matrices assambled, with
%              consistent units.
%   - time: Time vector, in [sec].
%   - ug_up: Updated vector with the acceleration time history, in [g].
%
%   u, v, a_t and Fs, have N rows, with N being equal to the number of
%   stories (N = length(M) = length(K)), and the same number of columns
%   than ug, representing the time variation of the response.
%
%   All the units must be consistent between the input parameters. The
%   output parameters will be also consistent.
%   
% Pablo Heresi
% 08-08-2015
%

N = length(wi);     % Number of stories
MaxIter = 20;

%% Obtain T and phi

% Note: first index corresponds to 1st floor, and last index to roof.
% M and K matrices
M = diag(wi)/g;         % Mass Matrix
K = ComputeK(k);        % Stiffness Matrix

% Eigenvalue analysis
[phi,w2] = eig(K,M);
w = sqrt(diag(w2));     % Undamped frequencies
[w,index] = sort(w);
T = 2*pi./w;            % Undamped periods

% Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0
sphi = phi;
for i = 1:N
    sphi(:,i) = phi(:,index(i))/ phi(end,index(i));
end
phi = sphi;             % Normalized modal shapes

% C matrix
Mi = diag(phi'*M*phi);
if size(Xi) == size(Mi)
    Ci = 2*Mi.*w.*Xi;
else
    Ci = 2*Mi.*w.*Xi';
end
C = (phi')^(-1)*diag(Ci)*phi^(-1);


%% Check stability of the method

Np = length(ug);     % Length of the record
time = 0:dt:(dt*(Np-1));

% If the time step is too large, display a warning.
if dt > T(1)/30
    disp(['Warning: The time step used (dt) for the ground motion "'...
        Name '" is greater than T_1/30. This is not recommended for representing'...
        ' the response correctly.']);
    disp(['A new dt = T/30 = ' num2str(T(1)/30) ' sec is used. In GM: ' Name]);
    
    dt_ = dt/ceil(30*dt/T(1));
    time_ = 0:dt_:time(end);
    ug = interp1(time,ug,time_);
    time = time_;
    dt = dt_;
    Np = length(ug);
end

%% Initial Calculations

r = ones(N,1);          % Note that this assumes horizontal excitation
P = -M*r*ug*g;          % Equivalent external load

dy = Fy./k;             % Yielding displacement of each story
dc = dcdy.*dy;          % Capping displacement of each story
Fmax = (1-a_s).*Fy+a_s.*k.*dc;      % Positive Strength Limit
Fmin = -Fmax;                       % Negative Strength Limit
LimMax = zeros(N,1);
LimMin = zeros(N,1);

% Initialize vectors
fs_st = zeros(N,Np);        % Story restoring force
fs = fs_st;                 % Total floor restoring force
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)];
Kt = K;                     % Initial tangent stiffness
kt = k;


u = zeros(N,Np);            % Relative displacement time history
v = u;                      % Relative velocity time history
a = u;                      % Relative acceleration time history
u(:,1) = do;                % Initial Displacement
v(:,1) = Vo;                % Initial Velocity
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1));  % Initial Relative Acceleration

% Constants
a1 = 4/dt^2*M + 2/dt*C;
a2 = 4/dt*M + C;
dt2 = dt/10;
a1_2 = 4/dt2^2*M + 2/dt2*C;
a2_2 = 4/dt2*M + C;

R = zeros(N,Np);            % Unbalanced force history
It = zeros(1,Np);           % Number of iterations
It(1) = 1;

%% Calculation for each time step
kt_prev = k;
i = 1;
kk = 0;
while i < size(u,2);

    u(:,i+1) = u(:,i);

    fs_st(:,i+1) = fs_st(:,i);
    fs(:,i+1) = fs(:,i);
    
    p_ = P(:,i+1) + a1*u(:,i) + a2*v(:,i) + M*a(:,i);
    
    % Newton-Raphson iterations
    j = 0;
    R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1);
    
    while sum(abs(R(:,i+1)) > tol)   &&  j < MaxIter+1
        Kt_ = Kt + a1;
        du = Kt_\R(:,i+1);
        u(:,i+1) = u(:,i+1) + du;
        
        % State determination
        [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
            (h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
        fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];
        
        Kt = ComputeK(kt);

        R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1);      % Unbalanced force
        j = j+1;                                % Increase # of iterations

        if j == MaxIter
            %disp(['Warning: Reached ' num2str(MaxIter) ' iterations. Convergence was not'...
            %    ' achieved in point i = ' num2str(i)])
        end
    end
    if [kt == kt_prev  ;  j < MaxIter+1]
        kk = kk+1;
        % After iterations are completed, compute relative vel and acc:
        It(i+1) = j;    % Total number of iterations
        v(:,i+1) = 2/dt*(u(:,i+1)-u(:,i)) - v(:,i);   % Rel Velocity
        a(:,i+1) = M\(P(:,i+1)-C*v(:,i+1)-fs(:,i+1)); % Relative acceleration
        kt_prev = kt;
        i = i+1;
    else    % Change in any stiffness or convergence not achieved
        dug2 = (ug(i+1)-ug(i))/10;

        time_int = (time(i)+dt2):dt2:(time(i+1)-dt2);
        if dug2 == 0
            ug_int = ug(i)*ones(1,9);
        else
            ug_int = ug(i)+dug2*(1:9);
        end
        
        time = [time(:,1:i) time_int time(:,i+1:end)];
        ug = [ug(:,1:i) ug_int ug(:,i+1:end)];

        P = -M*r*ug*g;          % Equivalent external load

        u = [u(:,1:i) zeros(N,9) u(:,i+1:end)];

        fs_st = [fs_st(:,1:i) zeros(N,9) fs_st(:,i+1:end)];
        fs = [fs(:,1:i) zeros(N,9) fs(:,i+1:end)];
        v = [v(:,1:i) zeros(N,9) v(:,i+1:end)];
        a = [a(:,1:i) zeros(N,9) a(:,i+1:end)];
        R = [R(:,1:i) zeros(N,9) R(:,i+1:end)];
        It = [It(:,1:i) zeros(1,9) It(:,i+1:end)];
        for i2 = 1:10
        	u(:,i+i2) = u(:,i+i2-1);
            fs_st(:,i+i2) = fs_st(:,i+i2-1);
            fs(:,i+i2) = fs(:,i+i2-1);
             
            p_ = P(:,i+i2) + a1_2*u(:,i+i2-1) + a2_2*v(:,i+i2-1) + M*a(:,i+i2-1);
                
                % Newton-Raphson iterations
                j = 0;
                R(:,i+i2) = p_ - fs(:,i+i2) - a1_2*u(:,i+i2);
                
                while sum(abs(R(:,i+i2)) > tol)   &&  j < MaxIter+1
                    Kt_ = Kt + a1_2;
                    du = Kt_\R(:,i+i2);
                    u(:,i+i2) = u(:,i+i2) + du;
                    
                    % State determination
                    [fs_st(:,i+i2),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
                        (h,Pi,u(:,i+i2-2),u(:,i+i2-1),u(:,i+i2),fs_st(:,i+i2-1),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
                    fs(:,i+i2) = [-diff(fs_st(:,i+i2));fs_st(end,i+i2)];
                    
                    Kt = ComputeK(kt);
                    
                    R(:,i+i2) = p_ - fs(:,i+i2) - a1_2*u(:,i+i2);      % Unbalanced force
                    j = j+1;                                % Increase # of iterations
                    
                    if j == MaxIter
                        %disp(['Warning: Reached ' num2str(MaxIter) ' iterations. Convergence was not'...
                        %    ' achieved in point i = ' num2str(i)])
                    end
                end
                % After iterations are completed, compute relative vel and acc:
                It(i+i2) = j;    % Total number of iterations
                v(:,i+i2) = 2/dt2*(u(:,i+i2)-u(:,i+i2-1)) - v(:,i+i2-1);   % Rel Velocity
                a(:,i+i2) = M\(P(:,i+i2)-C*v(:,i+i2)-fs(:,i+i2)); % Relative acceleration
        end
        i = i+10;
        kt_prev = kt;

    end

end
a_t = a/g + r*ug; 	% Absolute acceleration, in [g]

ug_up = ug;

end


function [K] = ComputeK(k)

if length(k) > 1
    k_aux = k(2:end);
    k_aux(end+1,1) = 0;
    K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);
else
    K = k;
end
end

    
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck)

du0 = [u0(1);diff(u0)];
du1 = [u1(1);diff(u1)];
du2 = [u2(1);diff(u2)];

Fmax2 = Fmax;
Fmin2 = Fmin;
LimMax2 = zeros(length(du1),1);
LimMin2 = zeros(length(du1),1);

fs = zeros(size(fs1));
fs1 = fs1 + Pi./h.*du1;
kt = k;
for i = 1:length(du1)
    fs(i) = fs1(i) + k(i)*(du2(i)-du1(i));
    
    if du2(i) > dc(i)
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)-dc(i));
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);
        ktEnv = a_c(i)*k(i);
    elseif du2(i) > -dc(i)
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);
        ktEnv = a_s(i)*k(i);
    else
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);
        NegLimEnv = (a_s(i)-1)*Fy(i)-a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)+dc(i));
        ktEnv = a_c(i)*k(i);
    end

    if fs(i) > min(Fmax(i),PosLimEnv)
        LimMax2(i) = 1;
        fs(i) = min(Fmax(i),PosLimEnv);
        if fs(i) == Fmax(i)
            kt(i) = 0;
        else
            kt(i) = ktEnv;
        end
    elseif fs(i) < max(Fmin(i),NegLimEnv)
        LimMin2(i) = 1;
        fs(i) = max(Fmin(i),NegLimEnv);
        if fs(i) == Fmin(i)
            kt(i) = 0;
        else
            kt(i) = ktEnv;
        end
    end
    
    if StrengthLimitCheck
        if du1(i) > dc(i) && du0(i) < du1(i) && du2(i) < du1(i) && LimMax(i)
            Fmax2(i) = fs1(i);
        end
        
        if du1(i) < -dc(i) && du0(i) > du1(i) && du2(i) > du1(i) && LimMin(i)
            Fmin2(i) = fs1(i);
        end
    end
    
    fs(i) = fs(i) - Pi(i)/h(i)*du2(i);
end
      
end                          
