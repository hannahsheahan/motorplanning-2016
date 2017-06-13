function  shadedTrajectory(x, y, colour, al)
% x: [x means; x standard errors]
% y: [y means; y standard errors]
% colour: colour of patch
% al: alpha transparency
%all numbers should be positively offset

%offset all mean numbers by 1000 just in case

x(1,:) = x(1,:) + 1000;
y(1,:) = y(1,:) + 1000;

%figure();

 for i=1:size(x, 2)
    
     xlow(i) = x(1,i) - x(2,i);
     xhigh(i) = x(1,i) + x(2,1);
     
     ylow(i) = y(1,i) - y(2,i);
     yhigh(i) = y(1,i) + y(2,i);
       
     % draw the patch on that datapoint (CCW from lower left corner)
     Xp = [xlow(i) xhigh(i) xhigh(i) xlow(i)];
     Yp = [ylow(i) ylow(i) yhigh(i) yhigh(i)];
     v = [xlow(i), ylow(i); xhigh(i), ylow(i); xhigh(i), yhigh(i); xlow(i), yhigh(i)];
     f = [1,2,3,4];
     patch('Faces',f,'Vertices', v, 'FaceColor', colour, 'EdgeColor', colour); hold on;
     %patch(Xp,Yp, 'FaceColor', colour, 'EdgeColor', 'none'); hold on;
     %alpha(al);
 end
 
   % Draw the mean line on top
   plot(x(1,:), y(1,:), 'k');
    
