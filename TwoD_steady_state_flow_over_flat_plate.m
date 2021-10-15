clear all
clc
nodes_x=11; %Number of nodes in the x direction
nodes_y=11; %Number of nodes in the y direction
x=linspace(0,1,nodes_x); %Grid
y=linspace(0,1,nodes_y);
L=1; %Length
no_of_iter=0; %Number of iterations
d_x=L/(nodes_x-1); %Distance between two grid points in x and y
d_y=L/(nodes_y-1);
beta=d_x/d_y; 
%Initializing the 'T' matrix that is the temperature distribution
T=ones(nodes_y,nodes_x); 
for i=1:nodes_x
    for j=1:nodes_y
        
        if(i==nodes_x)   %boundary condition for right boundary
            T(i,j)=0;
        end
        
        if(j==nodes_y)   %boundary condition for top boundary
            T(i,j)=0;
        end
    end
end
%'T' values initialized and plotting the countours
figure(1)
contourf(x,y,T);colormap(jet);
title(sprintf('Temperature Distribution at the beginning'));
%Fixing the temperature values in Told
Told=T;
%Enter the method Point Gauss-Siedal or Point successize over-relaxation
method=3;
%1-Jacobi
%2 - PGS
%3-PSOR
%Setting error to be a high value
error=9e9;
%Setting the tolerance value (This affects the accuracy and time for
%convergence)
tol=1e-06;
if(method==1)
    while(error>tol) %Outer loop checking for convergence
                for i=2:nodes_x-1
                    for j=2:nodes_y-1 % A function that goes along the x direction , row by row and applies the central differenced scheme to obtain the values of T(i,j)
                        T(i,j)=(1/(2+2*beta*beta))*(Told(i+1,j) + Told(i-1,j) + beta*beta*(Told(i,j+1) + Told(i,j-1))+(d_x*d_x));
                    end
                end
                
                for j=2:nodes_y-1 %Left Boundary
                    T(1,j)=(1/(2+2*beta*beta))*(Told(2,j) + Told(2,j) + beta*beta*(Told(1,j+1) + Told(1,j-1))+(d_x*d_x));
                end
                for i=2:nodes_x-1 %Bottom Boundary
                     T(i,1)=1/(2+2*beta*beta)*(Told(i+1,1) + Told(i-1,1) + beta*beta*(Told(i,2) + Told(i,2))+(d_x*d_x));
                end
                T(1,1)=1/(2+2*beta*beta)*(T(2,1) + T(2,1) + beta*beta*(T(1,2) + T(1,2))+(d_x*d_x)); %Corner Point
                error1=max(abs(Told-T)); %A matrix(1x6) of errors stored in error1
                error=max(error1);%Obtains the final value of the error.
                Told=T; %Fixing the Told values to the latest iteration
                no_of_iter=no_of_iter+1%Incrementing the number of iterations
                    


    end
end

if(method==2)
    while(error>tol) %Outer loop checking for convergence
                for i=2:nodes_x-1
                    for j=2:nodes_y-1 % A function that goes along the x direction , row by row and applies the central differenced scheme to obtain the values of T(i,j)
                        T(i,j)=(1/(2+2*beta*beta))*(Told(i+1,j) + T(i-1,j) + beta*beta*(Told(i,j+1) + T(i,j-1))+(d_x*d_x));
                    end
                end
                
                for j=2:nodes_y-1 %Left Boundary
                    T(1,j)=(1/(2+2*beta*beta))*(T(2,j) + T(2,j) + beta*beta*(Told(1,j+1) + T(1,j-1))+(d_x*d_x));
                end
                for i=2:nodes_x-1 %Bottom Boundary
                     T(i,1)=1/(2+2*beta*beta)*(Told(i+1,1) + T(i-1,1) + beta*beta*(T(i,2) + T(i,2))+(d_x*d_x));
                end
                T(1,1)=1/(2+2*beta*beta)*(T(2,1) + T(2,1) + beta*beta*(T(1,2) + T(1,2))+(d_x*d_x)); %Corner Point
                error1=max(abs(Told-T)); %A matrix(1x6) of errors stored in error1
                error=max(error1);%Obtains the final value of the error.
                Told=T; %Fixing the Told values to the latest iteration
                no_of_iter=no_of_iter+1%Incrementing the number of iterations
                    


    end
 end

 if(method==3)
    w=1.8; %Setting the value of 'w' 
    while(error>tol)%Outer loop checking for convergence
                for i=2:nodes_x-1 % A function that goes along the x direction , row by row and applies the central differenced scheme to obtain the values of T(i,j)
                    for j=2:nodes_y-1
                        T(i,j)= ((1-w)*Told(i,j))+(w/(2+2*beta*beta)*(Told(i+1,j) + T(i-1,j) + beta*beta*(Told(i,j+1) + T(i,j-1))+(d_x*d_x)));
                    end
                end
               
                for j=2:nodes_y-1 %Left Boundary
                     T(1,j)= ((1-w)*Told(1,j))+(w/(2+2*beta*beta)*(T(2,j) + T(2,j) + beta*beta*(Told(1,j+1) + T(1,j-1))+(d_x*d_x)));
                end

                for i=2:nodes_x-1 %Bottom Boundary
                    T(i,1)= ((1-w)*Told(i,1))+(w/(2+2*beta*beta)*(Told(i+1,1) + T(i-1,1) + beta*beta*(T(i,2) + T(i,2))+(d_x*d_x)));
                end
                 T(1,1)= ((1-w)*Told(1,1))+(w/(2+2*beta*beta)*(T(2,1) + T(2,1) + beta*beta*(T(1,2) + T(1,2))+(d_x*d_x))); %Corner Point
                  error1=max(abs(Told-T)); %A matrix(1x6) of errors stored in error1
                error=max(error1);%Obtains the final value of the error.
                Told=T;%Fixing the Told values to the latest iteration
                no_of_iter=no_of_iter+1 %Incrementing the number of iterations
    

    end
end
figure(2) %Plotting the final temperature countor.
contourf(x,y,T,10,'showText','on');colormap(jet);
colorbar;
title(sprintf('Temperature Distribution at the end'));