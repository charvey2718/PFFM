clear all
clc

%% Script parameters
global ls h b Gcint Gcbulk Gceff E nu alpha1 alpha2
E = 126e9; % Young's modulus (Pa)
nu = 0.3; % Poisson's ratio
ls = 20e-6; % length scale (m)
b = 0.05e-3; % half interface thickness (m)
Gcint = 281; % interface fracture toughness (N/m)
Gcbulk = 10*Gcint; % bulk fracture toughness (N/m)
h = b/12; % interface element size (m)
Gctol = 0.01; % Iterate until Gceff converges to within this tolerance
warning('off', 'all');

%% Yoshioka et al. (2021) as initial guess for ODE solver
Gceff = ((Gcbulk^2*h^2*exp((2*h)/ls) + Gcint^2*h^2*exp((2*h)/ls) + Gcbulk^2*h^2*exp((4*b)/ls) + Gcint^2*h^2*exp((4*b)/ls) + 4*Gcbulk^2*ls^2*exp((2*h)/ls) + 4*Gcint^2*ls^2*exp((2*h)/ls) + 4*Gcbulk^2*ls^2*exp((4*b)/ls) + 4*Gcint^2*ls^2*exp((4*b)/ls) - 2*Gcbulk^2*h^2*exp((h + 2*b)/ls) + 2*Gcint^2*h^2*exp((h + 2*b)/ls) + 8*Gcbulk^2*ls^2*exp((h + 2*b)/ls) + 8*Gcint^2*ls^2*exp((h + 2*b)/ls) - 2*Gcbulk*Gcint*h^2*exp((2*h)/ls) + 2*Gcbulk*Gcint*h^2*exp((4*b)/ls) + 8*Gcbulk*Gcint*ls^2*exp((2*h)/ls) + 8*Gcbulk*Gcint*ls^2*exp((4*b)/ls) - 4*Gcbulk^2*h*ls*exp((2*h)/ls) + 4*Gcint^2*h*ls*exp((2*h)/ls) + 4*Gcbulk^2*h*ls*exp((4*b)/ls) + 4*Gcint^2*h*ls*exp((4*b)/ls) - 48*Gcbulk*Gcint*ls^2*exp((h + 2*b)/ls) + 8*Gcint^2*h*ls*exp((h + 2*b)/ls) + 8*Gcbulk*Gcint*h*ls*exp((4*b)/ls) - 24*Gcbulk*Gcint*h*ls*exp((h + 2*b)/ls))^(1/2) + Gcbulk*h*exp(h/ls) + Gcint*h*exp(h/ls) - Gcbulk*h*exp((2*b)/ls) + Gcint*h*exp((2*b)/ls) - 2*Gcbulk*ls*exp(h/ls) + 2*Gcint*ls*exp(h/ls) - 2*Gcbulk*ls*exp((2*b)/ls) + 2*Gcint*ls*exp((2*b)/ls))/(2*h*exp(h/ls) + 2*h*exp((2*b)/ls) - 4*ls*exp(h/ls) + 4*ls*exp((2*b)/ls));
alpha1 = (exp((b - h/2)/ls)*(Gcbulk + Gceff))/(Gcbulk*exp((b - h/2)/ls) - Gcbulk*exp(-(b - h/2)/ls) + Gceff*exp((b - h/2)/ls) + Gceff*exp(-(b - h/2)/ls));
alpha2 =(2*Gceff*exp((b - h/2)/ls))/(Gcbulk*exp((b - h/2)/ls) - Gcbulk*exp(-(b - h/2)/ls) + Gceff*exp((b - h/2)/ls) + Gceff*exp(-(b - h/2)/ls));
fprintf('Yoshioka et al. (2021):\t%.2f\n', Gceff);

%% Method 1
[x_method1, s_method1, S_method1] = iterativelySolve(@f1, @bc, @guess_left, @guess_right, b, h, ls, Gcint, Gcbulk, Gctol);
fprintf('Method 1:\t%.2f\n', Gceff)
%% Method 2
[x_method2, s_method2, S_method2] = iterativelySolve(@f2, @bc, @guess_left, @guess_right, b, h, ls, Gcint, Gcbulk, Gctol);
fprintf('Method 2:\t%.2f\n', Gceff)
%% Method 3
[x_method3, s_method3, S_method3] = iterativelySolve(@f3, @bc, @guess_left, @guess_right, b, h, ls, Gcint, Gcbulk, Gctol);
fprintf('Method 3:\t%.2f\n', Gceff)

%% Functions
function [x, s, S] = iterativelySolve(f, bc, guess_left, guess_right, b, h, ls, Gcint, Gcbulk, Gctol)
    % The function f depends on Gceff, which is also being solved for. This
    % function therefore solves f iteratively, updating Gceff each time,
    % until Gceff is within tolerance.

    global Gceff

    xmesh_left = [linspace(0, b, 20000), linspace(b, 40*b, 20000)]; % starting from 0
    solinit_left = bvpinit(xmesh_left, guess_left);
    
    xmesh_right = [linspace(h/2, b, 20000), linspace(b, 40*b, 20000)];
    solinit_right = bvpinit(xmesh_right, guess_right);

    options = bvpset('RelTol', 1e-9, 'AbsTol', 1e-9);
    while true % iterate until Gceff is within tolerance
        sol_left = bvp4c(f, bc, solinit_left, options);
        sol_right = bvp4c(f, bc, solinit_right, options);
        x = [-flipud(sol_left.x'); sol_right.x'];
        s = [flipud(sol_left.y(1,:)'); sol_right.y(1,:)'];
        DeltaGceff = Gceff_from_sProfile(x, s, b, ls, Gcint, Gcbulk) - Gceff; % change in Gceff
        Gceff = Gceff + DeltaGceff; % update Gceff
        if abs(DeltaGceff) < Gctol % is Gceff within tolerance?
            break;
        end
    end

    dsdx = gradient(s)./gradient(x);
    S = 1/2*((1 - s).^2/ls + ls*dsdx.^2); % surface energy density function
end

function dsdx = f1(x, s, region)
    % The ODE to be solved from method 1.

    global ls Gcbulk Gceff
    dsdx(1, 1) = s(2);

    switch region
        case 1
            dsdx(2, 1) = (s(1).*(2*exp(-2*x/ls) + 1) - 1)/ls^2;
        case 2
            dsdx(2, 1) = (s(1).*(2*Gceff/Gcbulk*exp(-2*x/ls) + 1) - 1)/ls^2;
    end
end

function dsdx = f2(~, s, ~)
    % The ODE to be solved from method 2.

    global ls
    dsdx(1, 1) = s(2);
    k = ((1 - s(1)).^2/ls^2 + s(2).^2);
    dsdx(2, 1) = s(1)*(k + 1/ls^2) - 1/ls^2;
end

function dsdx = f3(x, s, region)
    % The ODE to be solved from method 3.

    global ls Gceff Gcint Gcbulk nu h
    
    % Calculate k from sharp crack theory
    switch region
        case 1
            k = (11 - 5*nu)/(16*ls*pi*x*s(1))*Gcint/Gceff;
        case 2
            k = (11 - 5*nu)/(16*ls*pi*x*s(1))*Gcint/Gcbulk;
    end

    % Limit k to kmax
    kmax = Gcint*(11 - 5*nu)/16/Gceff/ls/pi/(h/2);
    if k > kmax
        k = kmax;
    end
    
    dsdx(1, 1) = s(2);
    dsdx(2, 1) = s(1)*(k + 1/ls^2) - 1/ls^2;
end

function res = bc(YL, YR)
    % YL(i, j) is value of s or s' at left boundary of region j.
    % YR(i, j) is value of s or s' at right boundary of region j.
    % i = 1 for s, or i = 2 for s'.
    % For example, YR(2, 1) is s' on right of region 1.
    
    global Gcbulk Gceff
    res = [YL(1, 1); % s(0) = 0
        YR(1, 1) - YL(1, 2); % s1(b) = s2(b)
        Gceff*YR(2, 1) - Gcbulk*YL(2, 2); % first Weierstrass-Erdmann corner condition
        YR(1, 2) - 1]; % s(Inf) = 1;
end

function g = guess_left(x, ~)
    % Use Yoshioka et al.'s (2021) s profile without the element-size
    % effect as the initial guess for the left side.

    global alpha1 alpha2 ls
    if abs(x) > 0
        g = [1 - alpha2*exp(-(abs(x))/ls);
            (alpha2*exp((-abs(x))/ls))/ls];
    else
        g = [1 - alpha1*exp(-(abs(x))/ls) - (1 - alpha1)*exp((abs(x))/ls);
            (alpha1*exp((-abs(x))/ls))/ls + (exp(-(-abs(x))/ls)*(alpha1 - 1))/ls];
    end
end

function g = guess_right(x, ~)
    % Use Yoshioka et al.'s (2021) s profile with the element-size effect
    % as the initial guess for the right side.

    global alpha1 alpha2 h ls
    if abs(x) < h/2
        g = [0; 0];
    elseif abs(x) > h/2
        g = [1 - alpha2*exp(-(abs(x) - h/2)/ls);
            (alpha2*exp((h/2 - abs(x))/ls))/ls];
    else
        g = [1 - alpha1*exp(-(abs(x) - h/2)/ls) - (1 - alpha1)*exp((abs(x) - h/2)/ls);
            (alpha1*exp((h/2 - abs(x))/ls))/ls + (exp(-(h/2 - abs(x))/ls)*(alpha1 - 1))/ls];
    end
end

function Gceff = Gceff_from_sProfile(x, s, b, ls, Gcint, Gcbulk)
    % Apply the energy balance equation to compute Gceff.

    B = unique([x s],'rows');
    x = B(:,1); s = B(:,2);

    bulk_left  = x < -b;
    intf  = abs(x) <= b;
    bulk_right = x > b;
    
    ds_dx = gradient(s)./gradient(x);
    S = 1/2*((1 - s).^2/ls + ls*ds_dx.^2);

    Gceff = (Gcint - Gcbulk*trapz(x(bulk_left), S(bulk_left)) - Gcbulk*trapz(x(bulk_right), S(bulk_right)))/trapz(x(intf), S(intf));
end
