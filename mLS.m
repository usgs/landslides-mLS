% This script is provided as a supplementary material of a paper
% published in Earth Surface Processes and Landforms. The details of the 
% method followed in the given script is described in the corresponding paper.
% If you publish use this script or a its modified version please cite the 
% following paper:
% The entire link will be given when the paper will be ready to publish.

% The Pupose Of The Script And The Input Parametes 
% This script is provided for the accurate estimation of landslide-event
% magnitude. We used Matlab R2015b to test the given script. It requires 
% three input parameters: cutoff (smallest area that follows power law) 
% and beta (power-law exponent) values of Frequency-size distribution 
% of landslides, and a horizontal array (Area) with landslide sizes.  

% The power-law distribution can be captured in both cumulative and 
% non-cumulative FADs, and the power-law exponent for a non-cumulative 
% FAD (β) can be transferred to its cumulative equivalent, α, using 
% the relation α=β-1 (Guzzetti et al., 2002). In this code, the calculations 
% are carried out on the non-cumulative FAD.

% The cutoff is the landslide size where frequency-size distribution curve
% diverges form the power-law. In this code its unit is meter square.

% Beta value resfers to the slope of the frequency-size distribution. 
% It also called as power-law exponent (scaling parameter, β). 
% For most landslide inventories, non-cumulative power-law exponents occur
% in the range of 1.4–3.4, with a central tendency of 2.3–2.5 
% (Stark and Guzzetti, 2009; Van Den Eeckhaut et al., 2007).

% The unit of landslide sizes are in meter square.

% To obtain the cutoff and beta values the method suggested by Clauset 
% et al.(2009) should be used. The original script of Clauset et al. (2009)
% can be downloaded from the following link to calculate these parameters: 
% http://www.santafe.edu/~aaronc/powerlaws/ 

% The Output Of The Script
% When you run the mLS function, for the given sample data, the corresponding 
% mLS (landslide-event magnitude) will be obtained. As an output of this code, 
% a plot showing the frequency-area distribution of the given landslides and 
% the corresponding power-law fit are also obtained. mLS does not have any unit.

function [mLS]=mLS(Area,cutoff,beta)

% In the following lines, the bins are defined with an array. We took 2 as 
% the minimum bin size and we used increasing bin sizes. We increase the 
% bin widths while the landslide size increases, so that bin widths become
% approximately equal in logarithmic coordinates. To create a long array for 
% the bins we tentatively pick 120 for the size of x1 vector defined below. 
x1(1,1)=2;
for i=2:120
    x1(1,i)=x1(1,i-1)*1.2;
end
x2=log10(x1);

Freq=histc(Area,x1); %Frequency values are calculated for each bin 
s=size(x1);
s=s(1,2);
internal=zeros(1,s);

for i=2:s
     internal(1,i)=x1(1,i)-x1(1,i-1);
end
internal(1,1)=min(x1);
FD=Freq./internal;

x1_rev = abs(x1-cutoff);    % the index of value that is closest to cutoff value is identified along the x1 array 
[indexMidpoint indexMidpoint] = min(x1_rev);

x=x1(indexMidpoint:end); % the x (size bines) array for the frequeny-size distribution is defined 
y=FD(indexMidpoint:end); % the y (frequency densities) array for the frequeny-size distribution is defined 

if beta>0           % beta value have to be negative
    beta=-1*beta;
end

constant=y(1,1)/cutoff^beta;     % The c constant is calculated along the power-low where x=cutoff
fit_y=constant*x1.^beta;        % Frequency-density values calculated for the defined power-law fit

midx=10^((log10(max(Area))+(log10(cutoff)))/2);     % The x and y values at mid-point location is read along the power-law fit
midy=constant*midx^beta;

Refmidx=2.503413777225013e+04;    % X value for the mid point of the Haiti (reference point) inventory
Refmidy=0.005717876265502;        % Y value for the mid point of the Haiti (reference point) inventory
ac=Refmidy/(23559*Refmidx^beta);  % he c' constant (as) is calculated here for the mid-point of 
                                  % the Haiti inventory as a reference point where mLS=log(23559)=4.37 
                             
mLS=log10((midy/(ac*midx^(beta)))); % mLS is calculated in this line

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
