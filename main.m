%% MAIN FUNCTION

function main()

% User parameters
plot_on = 1;
max_iterations = 1000;

L = 1.0;
H = 0.5;
Nx = 3;
Ny = 3;



% Define system

% Make initial guess
Temperature_vector = rand(Nx*Ny,1);  

% Initialize grid
[gridx,gridy]=generate_grid(L,H,Nx,Ny);

delta(3,3,'N','P',gridx,gridy)

% Build load vector
load_vector=zeros(Nx*Ny,1);

% Build mass matrix
Mass_matrix = zeros(Nx*Ny);

for i = 1:Nx*Ny
    % Setting mass matrix
    Mass_matrix(i,i)=a_P(i);
    if i <= Nx*(Ny-1) %Not North edge 
        Mass_matrix(i,i+Nx)=a_N(i);
    end

    if mod(i,Nx)~= 0 %Not East edge
        Mass_matrix(i,i+1)=a_E(i);
    end

    if mod(i,Nx)~= 1 %Not West edge
        Mass_matrix(i,i-1)=a_W(i);
    end        

    if i > Nx %Not South edge 
        Mass_matrix(i,i-Nx)=a_S();
    end           
end
    

% Solve system
%    Temperature_vector=Mass_matrix\Source_vector; % This is not allowed...
for step = 1:max_iterations;
    Temperature_vector=Gauss_Siedel_Step(Mass_matrix,load_vector,Temperature_vector);
end
% Visualize system
if plot_on == 1
    subplot(1,2,1)
    imagesc(Mass_matrix);
    subplot(1,2,2)
    plot_field_1d(Temperature_vector,Nx,Ny)
end


end


%% Grid functions
function [gridx,gridy] = generate_grid(L,H,Nx,Ny)
%[gridx,gridy] = generate_grid(L,H)
    % GENERATE A GRID
    % TODO:
    % This is probably VERY buggy!
    % ?Make nonlinear grid
    
    cell_lengthx = L/Nx;
    cell_lengthy = H/Ny;
    
    
    %Center positions
    centerx=linspace(cell_lengthx/2,L-cell_lengthx/2,Nx)';
    centery=linspace(cell_lengthy/2,H-cell_lengthy/2,Ny)';
    % Cell sizes
    lengthx=ones(Nx,1)*cell_lengthx;
    lengthy=ones(Ny,1)*cell_lengthy;
    % Final grids
    gridx=[centerx,lengthx];
    gridy=[centery,lengthy];
end

function [x,y]= position(cell_x,cell_y,compass,gridx,gridy)
% [x,y]=POSITION(cell,compass,gridx,gridy)
    % Function to get position of point
    % Working
    %
    % x         float               Position x in grid 
    % y         float               Position y
    %
    % cell      integer             Cell index
    % compass   P,N,E,S,W,n,e,s,w   Point in cell
    % gridx     [Nx*2]              [cell1_center_x, cell1_length_x; cell2_center_x ...]
    % gridy     -||-                [cell1_center_y, cell1_length_y; cell2_center_y ...]
    
    assert(ischar(compass),'compass is not a string')
    
    x = gridx(cell_x,1);
    y = gridy(cell_y,1);    
    if strcmp(compass, 'P')
        
    elseif strcmp(compass, 'n')
        y = y+gridy(cell_y,2)/2;
    elseif strcmp(compass, 'e')
        x = x+gridx(cell_x,2)/2;
    elseif strcmp(compass, 's')
        y = y-gridy(cell_y,2)/2;
    elseif strcmp(compass, 'w')
        x = x-gridx(cell_x,2)/2;
    elseif strcmp(compass, 'N')
        y = y+gridy(cell_y,2);
    elseif strcmp(compass, 'E')
        x = x+gridx(cell_x,2);
    elseif strcmp(compass, 'S')
        y = y-gridy(cell_y,2);
    elseif strcmp(compass, 'W')
        x = x-gridx(cell_x,2);
    else 
        x = NaN;
        y = NaN;
        error('Non valid compass direction')
    end
end

function [length]=delta(cell_x,cell_y,compass1,compass2,gridx,gridy)
    [x1,y1] = position(cell_x,cell_y,compass1,gridx,gridy);
    [x2,y2] = position(cell_x,cell_y,compass2,gridx,gridy);
    dx=x1-x2;
    dy=y1-y2;
    length=sqrt(dx.^2+dy.^2);
end

%% Assembly functions

function [i]=index_1d(x,y,Nx,~)
% [i]=index_1d(x,y,Nx,~)
    % Get the 1D-index from the x- and y-index
%     assert(isinteger(x))
%     assert(isinteger(y))
%     assert(isinteger(Nx))
    i=(y-1)*Nx+x;
end

function [x,y]=index_2d(i,Nx,~)
% [x,y]=index_2d(i,Nx,~)
    % Get the 2D-index from the i-index.
%     assert(isinteger(i))
%     assert(isinteger(Nx))
%     
    x=mod(i,Nx);
    if x==0
        x=Nx;
    end
    y=(i-x)/Nx+1;
end

function [a]=a_P(i)

a=1;
end

function [a]=a_N(i)

a=1;
end

function [a]=a_E(i)

a=1;
end

function [a]=a_S(i)

a=1;
end

function [a]=a_W(i)

a=1;
end

function [b]=b_bulk(i)

b=1;
end

function [b]=b_N(i)

b=1;
end

function [b]=b_E(i)

b=1;
end

function [b]=b_S(i)

b=1;
end

function [b]=b_W(i)

b=1;
end

%% Gauss-Siedel solver

function [T]=Gauss_Siedel_Step(A,b,T)
    assert(length(b)==length(T))
    N=length(b);
    for i = 1:N
        sigma=0;
        for j= 1:N
            if j~=i
                sigma = sigma + A(i,j).*T(j);
            end
        end
        T(i)=1./A(i,i).*(b(i)-sigma);
    end
end

%% Plotting functions

function [map]=field_2_2d(field,Nx,Ny)
    assert(length(field)==Nx*Ny);
    map=zeros(Nx,Ny);
    for i = 1:length(field)
       [x,y]=index_2d(i,Nx);
       map(x,y)=field(i);
    end
end

function plot_field_1d(field,Nx,Ny)
    map = field_2_2d(field,Nx,Ny);
    imagesc(map)
end
