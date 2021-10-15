%%Consider a hot plate of finite thickness 2L, which is suddenly exposed to
%%a fluid at T_infinity. Initial temperature of the plate is T_i. Heat
%%transfer coefficient is large. Find the temperature distribution of the
%%plate as a function of space and time.%%
clc
clear;
alpha = 1; % k/(rho*c)
L = 1;  % length of domain
T = 1; % maximum time upto which the simulation needs to be run


Nt = 400;   %no. of time-step to reach max time
Dt = T/Nt;  % delta t = time step size
Nx=10;      % No. grid point-1
Dx=L/Nx;    % grid size

b=alpha*Dt/(Dx*Dx); % beta = alpha*delat t/(delat x)^2

disp(b)


for i = 1:Nx+1      % Loop for Intial Condition
    x(i)=(i-1)*Dx;  % progression in space from i=1 to Nx+1
    th(i,1)=1;      % value of all theta at initial time (t=1)
end

for k=1:Nt+1        % Loop for Boundary Condition
    th(1,k)=0;      % value @ left boundray i.e theta @ i=1,at all time
    t(k)=(k-1)*Dt;  % progression in time from time t=1 to Nt+1
end

for k=1:Nt          % time loop
    for i=2:Nx      % spatial loop
        th(i,k+1)=th(i,k)+ (b*(th(i-1,k)+th(i+1,k)-(2*th(i,k))));       % explicit equation for internal nodes
    end
    th(Nx+1,k+1)= th(Nx+1,k)+ (b*(2*th(Nx,k) - 2*th(Nx+1,k)));          % Explicit equation for right boundary node
end

figure (1)
plot(x,th(:,1),'-b',x,th(:,Nt/20),'--g',x,th(:,Nt/40),'-r',x,th(:,Nt),':b')


