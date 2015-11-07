% MAIN FUNCTION

function main()
% User parameters

L = 1.0;
H = 0.5;
Nx = 10;
Ny = 10;

% Define system
    % Initialize grid
    [gridx,gridy]=generate_grid(L,H,Nx,Ny);
    
    delta(3,3,'N','P',gridx,gridy)
    % Build mass matrix
    Mass_matrix = rand(Nx*Ny);
    
    % Make initial guess
    Temperature_vector = ones(Nx*Ny,1);
% Solve system

% Visualize system

contourf(Mass_matrix);
disp(Temperature_vector');
disp(gridx')

end


% Grid functions
function [gridx,gridy] = generate_grid(L,H,Nx,Ny)
%[gridx,gridy] = generate_grid(L,H)
    % GENERATE A GRID
    % TODO:
    % This is probably VERY buggy!
    % ?Verify that this is a good choice of grid (add L/(2Nx)?)
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


