function [rho1, theta1, rho2, theta2] = detect_lanes (~)

setN = [1 2]; % Number of sets in train
totalFramesInSets = [250 232]; % Total number of frames in set1 and set2

%reads and converts image to grayscale
for i = setN  
    for j = 1:totalFramesInSets(i)
        % load each frame in the set
        image_name = sprintf('set%d_f%05d.png',i,j-1);
        img = imread(image_name);

%converts image to grayscale
%img=int8(img);
img=rgb2gray(img);
%img = wiener2(img); %removes noise from image
img = imgaussfilt(img,6); %guassian filter
%img=im2bw(img);
img_edge=edge(img,'canny');

tFreq=0.01; %theta frequency
% tFreq=1/2000; %theta frequency
% 
[w,h]=size(img_edge);

rhoLimit=norm([w,h]);

% %first lane
% %theta1 = pi/4; %2.5185;
% theta1=rhoLimit;
% rho1=100;
% %rho1 = x*cos(theta1) + y*sin(theta1);
% 
% %second lane
% %theta2 = -pi/4; %-0.6646;
% theta2=-rhoLimit;
% rho2=200;
% %rho2 = x*cos(theta2) + y*sin(theta2);
% 
rho=(-rhoLimit:1:rhoLimit);
% rho=(theta2:1:theta1); 
theta=(0:tFreq:pi);

numThetas = numel(theta); %number of arrays
houghSpace = zeros(numel(rho),numThetas);

%hough transform
% for wi = 1:w
%     for hj = 1:h
%         if image(hj, wi)==1 
%             for theta_index=1:numThetas
%                 th = theta(theta_index);
%                 r  = wi * cos(th) + hj * sin(th);
%                 rho_index = round(r + num_rhos/2);                      
%             houghSpace(rho_index, theta_index)+1;
%             end
%         end
%     end
% end

accum_x=ceil((theta(3)-theta(1))/theta(2))+1;

% rho_max=ceil(norm(size(img_edge)));
% rho_min=0;
% d_rho=rho_max/accum_x;
% rho=[rho_min, d_rho, rho_max];
accum_y=accum_x;

accum=zeros(accum_x,accum_y);

%find peaks in the hough transform
r=[];
c=[];
[maxCol,rowNum]=max(houghSpace);
% [rows,colms]=size(img_edge);
difference = 25;
thresh = max(max(houghSpace))-difference;
for m = 1:size(maxCol, 2)
   if maxCol(m)>thresh
       c(end+1)=m;
       r(end+1)=rowNum(m);
   end
end

% houghPeak=plot(theta(c),rho(r),'rx'); %plots the peaks on hough transformed image

for n = 1:size(c,2)
    th = theta(c(n));
    rh = rho(r(n));
    m = -(cos(th)/sin(th));
    b = rh/sin(th);
    x = 1:maxCol;
    plot(x, m*x+b);
end

% lines = houghlines(img_edge,theta,rho,houghPeak,'FillGap',5,'MinLength',30);
%     
% %     figure, imshow(img), hold on
%         max_len = 0;
%         for k = 1:length(lines)
%            xy = [lines(k).point1; lines(k).point2];
%            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%            % Plot beginnings and ends of lines
%            plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%            plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% 
%            % Determine the endpoints of the longest line segment
%            len = norm(lines(k).point1 - lines(k).point2);
%            if ( len > max_len)
%               max_len = len;
% %               xy_long = xy;
%            end
%         end

% for x=1:size(img_edge,1)
%   for y=1:size(img_edge,2)
%       theta_i = 1;
%       for t=theta(1):theta(2):theta(3)
% %         r = x * cos(t) + y * sin(t);
%         r = cos(t) + sin(t);
%         rho_i = ceil((r-rho(1))/rho(2));
%         accum(theta_i,rho_i) = accum(theta_i,rho_i) + 1;
%         theta_i = theta_i + 1;
%       end
%   end
% end
% 
% theta1=theta;
% rho1=rho_i;
% theta2=-theta;
% rho2=rho_i/2;

%finds edge pixels
[xI,yI] = find(img_edge);

%accumulator array
edgePix=numel(xI);
acc= zeros(edgePix,numThetas);

cosine = (0:w-1)'*cos(theta);
sine = (0:h-1)'*sin(theta);

acc((1:edgePix),:)=cosine(xI,:)+sine(yI,:);

max_val = max(acc(:));

thold = .01*max_val;

loc = acc >= thold;

count=sum(loc(:));

params = zeros(count,2);
index = 1;

for x=1:size(loc,1)
  for y=1:size(loc,2)
    if loc(x,y) > 0
      t = theta(1)+(x-1)*theta(2);
      r = rho(1)+(y-1)*rho(2);
      params(index,:) = [t,r];
      index = index + 1;
      
      theta1=t;
      rho1=r;
%       rho1=pi/4;
%       rho2=rho1/2; % 139;
      rho2=rho1/2;
%       theta2 = -pi/4; %-0.6646;
      theta2=-t;
    end
  end
end

imshow(img);

%votes
    for k = (1:numThetas)
            houghSpace(:,k) = hist(acc(:,k),rho);
    end
    end
end
