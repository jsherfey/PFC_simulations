function varargout = plot_CI(x,y,err,c,marker)
% purpose: Linear plot with continuous confidence/error boundaries.
%
%   plot_CI(X,Y,ERR,C) plots the graph of vector X vs. vector Y with
%   'continuous' confidence/error boundaries specified by the vector
%   ERR in the color specified by C.
%
%   H = plot_CI(...) returns a vector of line handles.
%
%   Example:
%      x = 1:0.1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      plot_CI(x,y,e,[1 0 0])
%   draws symmetric continuous confidence/error boundaries of unit standard deviation.
%
if nargin<5, marker='o'; end
if nargin<4, c='k'; end
if nargin<3, err=0; end
if (nargin<2)
    disp('ERROR: not enough input arguments!');
    return;
end

% transpose vectors (x,y,err) if incorrect dimensions
if size(x,1)>size(x,2), x=x'; end
if size(y,1)>size(y,2), y=y'; end
if size(err,1)>size(err,2), err=err'; end

z1 = y + err;
z2 = y - err;

p = plot(x,y,x,z1,x,z2);    YLIM = get(gca,'YLim');    delete(p);

X=[x fliplr(x)]';Y=[z1 fliplr(z2)]';
pa1=patch(X,Y,c,'FaceAlpha',0.4,'LineStyle','none');set(gca,'XLim',[min(x) max(x)]);
set(get(get(pa1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
hold on;
p = plot(x,y,[c marker '-'],'LineWidth',2);
hold off;

set(gca,'Layer','top');

H = [p, pa1];

if (nargout>0) varargout{1} = H; end;
