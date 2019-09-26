function contour_fill(x,y,z,v,c,flag)
%   function contour_fill(x,y,z,v,c,[flag])
% Fill the contour region z>v with a specified color. Similar to
% contourf(x,y,z,[v v]), except that the contour region is filled with
% the color c rather than using the current colormap.
% This script is a bit finicky and does not always work. Try using flag=1
% if it fails with default (flag=0).
% Example:
%   z = peaks; [x,y]=meshgrid(1:size(z,2),1:size(z,1));
%   clf, contour_fill(x,y,z,1,[0 0 1])
%   hold on, contour_fill(x,y,z,3,[0 0.5 0])
%   contour_fill(x,y,z,5,[1 0 0]),
%   contour(x,y,z,[-5:5]), hold off
% Ian Eisenman, 2008

if nargin~=5 && nargin~=6, disp('contour_fill(x,y,z,v,c)'), return; end

if nargin==5, flag=0; end

if max(size(x))>1 && diff(x(1:2))~=0 % x,y are arrays and transposed from expected orientation
    x=x'; y=y'; z=z';
end
if max(size(x))>1 && diff(x(1:2))~=0 % not a uniform grid
    disp('ERROR: Not a uniform grid')
    return
end

hold_status=ishold;
hold on

if flag==0 % one method that often works: switch "if 1" to "if 0" to use other method
    a=contourc(x(1,:),y(:,1),z,[v v]); % compute contour lines
    ki=1;
    while length(a)>1
        k=(ki+1):(ki+a(2,ki));
        %plot(a(1,k),a(2,k),'k')
        fill(a(1,k),a(2,k),c,'edgecolor','none')
        ki=ki+a(2,ki)+1;
        if ki>length(a), a=0; end % stop when lines have all been drawn
    end
else % another method that occasionally works
    a=contourc(x(1,:),y(:,1),z,[v v]); % compute contour lines
    ki=1; X=[]; Y=[];
    while length(a)>1
        k=(ki+1):(ki+a(2,ki));
        X=[X a(1,k)]; Y=[Y a(2,k)];
        if ki==1, X=fliplr(X); Y=fliplr(Y); end
        ki=ki+a(2,ki)+1;
        if ki>length(a), a=0; end % stop when lines have all been drawn
    end
    fill(X,Y,c,'edgecolor','none')
end

hold off
if hold_status, hold on, end
