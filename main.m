% MAIN FUNCTION

function main()

% Define system
    % Initialize grid
    
    % Build mass matrix
    
    % Make initial guess

% Solve system

% Visualize system

end


% Grid functions
function [gridx,gridy] = generate_grid(L,H)

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




