function [rho1, theta1, rho2, theta2] = detect_lanes (img)

% setN = [1 2]; % Number of sets in train
% totalFramesInSets = [250 232]; % Total number of frames in set1 and set2

%reads and converts image to grayscale
% for i = setN
%     for j = 1:totalFramesInSets(i)
    % load each frame in the set
%     image_name = sprintf('set%d_f%05d.png',i,j-1);
%     image_name = sprintf('set1_f00000.png',i,j-1);
    img_orig = imread('set1_f00000.png');
    
    img_gray=rgb2gray(img_orig);
    img = imgaussfilt(img_gray,2); %guassian filter
    img_edge=edge(img,'canny');
%     [r,c]=find(img_edge);
%     edgeIm=[c,r];
    
%     maxDist=1.2;
%     numTheta=200;
%     ns=max(edgeIm(:))*maxDist;
%     row=edgeIm(2,:);
%     col=edgeIm(1,:);
%     ts=[0:pi/numTheta:pi-pi/numTheta]; %range of line
%     
%     %cos ands in of angles
%     cTheta=cost(ts);
%     sTheta=sin(ts);
%     
%     s=row*cTheta'+col*sTheta';
%     
%     smin=min(s(:));
%     smax=max(s(:));
%     m=(ns-1)/(smax-smin);
%     b=(smax-ns*smin)/(smax-smin);
%     rs=round(m*s+b);
%     
%     h=[];
%     for k=1:ns
%         isEq=(rs==k);
%         h(:,k)=sum(isEq);
%     end
%     
    


    [H,theta,rho]=hough(img_edge);
    p=houghpeaks(H,10,'threshold',ceil(0.1*max(H(:))));
    lines = houghlines(img_edge,theta,rho,p,'FillGap',20,'MinLength',40);
    
    figure, imshow(img_orig), hold on
    max_len = 0;
        for k = 1:length(lines)
           xy = [lines(k).point1; lines(k).point2];
           plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

           % Plot beginnings and ends of lines
           plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
           plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

           % Determine the endpoints of the longest line segment
           len = norm(lines(k).point1 - lines(k).point2);
           if ( len > max_len)
              max_len = len;
              xy_long = xy;
           end
        end
    % highlight the longest line segment
    plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

%     rho1 = 300;
%     theta1 = theta; %2.5185;
%     rho2 = 200; % 139;
%     theta2 = -theta; %-0.6646;

end
