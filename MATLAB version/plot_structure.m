%% STRUCTURE FIGURE
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

% Input data:
% nodes: text file with information about structure nodes
    ... column 1 and 2: node X and Y coordinate
    ... column 3, 4 and 5: node boundary conditions (horizontal displacement, 
%             vertical displacement, rotation (if = 1 movement is restricted)
%             obs., truss elements have no rotation, so there is no columns 
%             5 and 8
% bars: text file with information about structure bars
   ... column 1 and 2: N1, N2 = element nodes


function plot_structure(nodes,bars,element_type)

    for i = 1:size(bars,1)
        N1 = bars(i,1);
        N2 = bars(i,2);

        x1 = nodes(N1,1);
        y1 = nodes(N1,2);
        x2 = nodes(N2,1);
        y2 = nodes(N2,2);
        
        x = [x1 x2];
        y = [y1 y2];
        
        xmax = max(nodes(:,1));
        if xmax == 0
            xmax = 1;
        end
        ymax = max(nodes(:,2));
        if ymax == 0
            ymax = 1;
        end
        
        plot(x,y,'black')
        hold on
        scatter(x,y,'filled','black')
        xlim([-(0.2*xmax) (xmax+(0.2*xmax))])
        ylim([-(0.1*ymax) (ymax+(0.1*ymax))])
        hold on

    end
    
    for i = 1:size(nodes,1)
       
        X = nodes(i,1);
        Y = nodes(i,2);
        RX = nodes(i,3);
        RY = nodes(i,4);
        
        if RX == 1   % zero X-displacement boundary conditions
           scatter(X,Y,220,'>','r') 
           hold on
        end
        
        if RY == 1   % zero Y-displacement boundary conditions
           scatter(X,Y,220,'^','r') 
           hold on
        end
        
        if strcmp(element_type,'timoshenko beam') || strcmp(element_type,'plane frame')
            RZ = nodes(i,5);
            if RZ == 1   % zero Z-rotation boundary conditions
                scatter(X,Y,500,'square','filled','r') 
                hold on
            end
        end
        
    end
   
end