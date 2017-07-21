% This script is provided for the accurate estimation of landslide-event
% magnitude. This script requires three input parameters: cutoff 
% (smallest area that follows power law) and beta (power-law exponent) 
% values of Frequency-size distribution of landslides, % and a horizontal 
% array (Area) with landslide sizes.   

% To obtain the cutoff and beta values the method suggested by Clauset 
% et al.(2009) should be used. The original script of Clauset et al. (2009)
% can be downloaded from the following link to calculate these parameters: 
% http://www.santafe.edu/~aaronc/powerlaws/ 

% When you run the mLS function, for the given sample data, the corresponding 
% mLS (landslide-event magnitude) will be obtained. As an output of this code, 
% a plot showing the frequency-area distribution of the given landslides and 
% the corresponding power-law fit are also obtained.

function [mLS]=mLS(Area,cutoff,beta)

% In the following lines, the bins are defined with an array. We took 2 as 
% the minimum bin size and we used increasing bin sizes. We increase the 
% bin widths while the landslide size increases, so that bin widths become
% approximately equal in logarithmic coordinates.
x1(1,1)=2;
for i=2:120
    x1(1,i)=x1(1,i-1)*1.2;
end
x2=log10(x1);

Freq=histc(Area,x1); %Frequency values are calculated for each bin 
s=size(x1);
s=s(1,2);
int=zeros(1,s);

for i=2:s
     int(1,i)=x1(1,i)-x1(1,i-1);
end
int(1,1)=min(x1);
FD=Freq./int;

tmpm1 = abs(x1-cutoff); % the index values for the mid-point is calculated for the inventory  
[idxm1 idxm1] = min(tmpm1);

x=x1(idxm1:end);
y=FD(idxm1:end);

if beta>0           % beta value have to be negative
    beta=-1*beta;
end

constant=y(1,1)/cutoff^beta;     % The c constant is calculated along the power-low where x=cutoff
fit_y=constant*x1.^beta;        % Frequency-density values calculated for the defined power-law fit

midx=10^((log10(max(Area))+(log10(cutoff)))/2);     % The x and y values at mid-point location is read along the power-law fit
midy=constant*midx^beta;

Nmidx=2.503413777225013e+04; % T
Nmidy=0.005717876265502;
as=Nmidy/(23559*Nmidx^beta); % he c' constant (as) is calculated here for by taking the mid-point of 
                             % the Haiti inventory as a reference point where mLS=log(23559)=4.37 
                             
mLS=log10((midy/(as*midx^(beta)))); % mLS is calculated in this line

% A graph showing the frequency-area distribution of the given landslides 
% and the corresponding power-law fit are plotted.

loglog(x1,fit_y,'-','LineWidth',2,'Color','r');hold on
loglog(x1,FD,'ok')
axis([1 1.E+7 1.E-6 1000])
set(get(gca,'Xlabel'),'string','Landslide Area (m^2)','FontSize',12, 'FontUnits','points','FontWeight','normal')
set(get(gca,'Ylabel'),'string','Frequency Density (m^-^2)','FontSize',12, 'FontUnits','points','FontWeight','normal')
str={['\beta = ',num2str(beta)];['mLS = ',num2str(mLS)]};
text(midx*2,midy*2,str)
end
