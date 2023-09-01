function [rho_intersect,e_intersect,hasintersect] = intersection(U_old, U_new, curve)

% U_old: U = [rho, rho*e, T] before intersection
% U_new: U = [rho, rho*e, T] before intersection
% curve: [rho, e] pairs, 

% line 1:
XY1 = [U_old(1) U_old(2)/U_old(1) U_new(1) U_new(2)/U_new(1)];
% curve made of line segments
% note, input is of the form [rhog eg]
XY2 = [curve(1:end-1,:) curve(2:end,:)];

%% possible plot:
plot(XY1([1 3]),XY1([2 4]),'s-')
hold on
plot(curve(:,1),curve(:,2))

%%

out = lineSegmentIntersect(XY1,XY2);
ind_intersect = find(out.intAdjacencyMatrix);
if (isempty(ind_intersect)) % no intersection, take new point
    hasintersect = false;
    rho_intersect = NaN;
    e_intersect = NaN;
else
    hasintersect = true;
    rho_intersect = out.intMatrixX(ind_intersect);
    e_intersect   = out.intMatrixY(ind_intersect);
end

