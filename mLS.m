% This script is provided as a supplementary material of a paper
% published in Earth Surface Processes and Landforms. The details of the 
% method followed in the given script is described in the corresponding 
% paper. If you publish use this script or a its modified version please 
% cite the following paper:
% Tanyas, H., K.E. Allstadt, and C.J. van Westen, 2018, 
% An updated method for estimating landslide-event magnitude, Earth Surface 
% Processes and Landforms. DOI: 10.1002/esp.4359

% The Pupose Of The Script And The Input Parametes 
% This script is provided for the accurate estimation of landslide-event
% magnitude. We used Matlab R2015b to test the given script. It basically
% requires three input parameters to estimate landslide-event magnitude: 
% cutoff smallest area that follows power law)and beta (power-law exponent)
% values of frequency-size distribution (FAD) of landslides, and a horizontal
% array (Area) with landslide sizes. This script can also calculate the
% uncertainty in landslide-event magnitude if uncertainties in cutoff and
% beta values are given as input parameters.

% The power-law distribution can be captured in both cumulative and 
% non-cumulative FADs, and the power-law exponent (beta) for a non-cumulative 
% FAD can be transferred to its cumulative equivalent, alpha, using the relation 
% alpha=beta-1 (Guzzetti et al., 2002). In this code, the calculations are carried 
% out on the non-cumulative FAD.

% The cutoff is the landslide size where frequency-size distribution curve
% diverges form the power-law. In this code its unit is meter square.

% Beta value resfers to the slope of the frequency-size distribution. It also
% called as power-law exponent (scaling parameter, beta). For most landslide
% inventories, non-cumulative power-law exponents occur in the range of 
% 1.4–3.4, with a central tendency of 2.3–2.5 (Stark and Guzzetti, 2009;
% Van Den Eeckhaut et al., 2007).

% The unit of landslide sizes are in meter square.

% To obtain the cutoff & beta values and thier uncertainties the method
% suggested by Clauset et al.(2009) can be used. The original scripts 
% (plfit.m and plvar.m) of Clauset et al. (2009) can be downloaded from the
% following link to calculate these parameters: 
% http://www.santafe.edu/~aaronc/powerlaws/ 

% The Output Of The Script
% When you run the mLS function, for the given sample data, the corresponding 
% mLS (landslide-event magnitude) will be obtained. If the uncertainties in 
% cutoff (cutoff_error) and beta (beta_error) values are provided, the
% uncertainty in mLS (error) is also estimated by the script. As an output
% of this code, a plot showing the frequency-area distribution of the given 
% landslides and the corresponding power-law fit are also obtained.
% mLS does not have any unit.

function [mLS,error]=mLS(Area,cutoff,beta,beta_error,cutoff_error)

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
beta_stored=beta;

constant=y(1,1)/cutoff^beta;     % The c constant is calculated along the power-low where x=cutoff
fit_y=constant*x1.^beta;        % Frequency-density values calculated for the defined power-law fit
fit_y_stored=fit_y;             

midx=10^((log10(max(Area))+(log10(cutoff)))/2);     % The x and y values at mid-point location is read along the power-law fit
midy=constant*midx^beta;

Refmidx=4.876599623713225e+04;    % X value for the mid point of the Northridge (reference point) inventory
Refmidy=8.364725347860417e-04;    % Y value for the mid point of the Northridge (reference point) inventory
ac=Refmidy/(11111*Refmidx^beta);  % the c' constant (as) is calculated here for the mid-point of 
                                  % the Northridge inventory as a reference point where mLS=log(11111) 
                             
mLS=log10((midy/(ac*midx^(beta)))); % mLS is calculated in this line
mLS_stored=mLS;

% Uncertainty in mLS will be calculated if the required inputs are given
% To do that Monte-Carlo simulation is run 10,000 times 

if exist('beta_error')==1 && exist('cutoff_error')==1  
    beta_interval_N=((beta+beta_error)-(beta-beta_error))/(499); %Number of elements to create an array including 500 uniformly distributed beta values is defined 
    beta_interval=(beta-beta_error):beta_interval_N:(beta+beta_error);  %An array including 500 uniformly distributed beta values is defined     

    cutoff_min=(cutoff-cutoff_error); 
        if cutoff_min<=0
             cutoff_min=2;
        else
             cutoff_min=cutoff_min;
        end
    cutoff_max=(cutoff+cutoff_error);
    cutoff_interval_N=((cutoff_max)-(cutoff_min))/(499);    %Number of elements to create an array including 500 uniformly distributed cutoff values is defined 
    cutoff_interval=(cutoff_min):cutoff_interval_N:(cutoff_max); %An array including 500 uniformly distributed cutoff values is defined   

    beta_mean=mean(beta_interval); %Mean of beta values is identified
    beta_std=std(beta_interval);    %Standard deviation of beta vaues is identified
    cutoff_mean=mean(cutoff_interval);  %Mean of cutoff values is identified
    cutoff_std=std(cutoff_interval);    %Standard deviation of cutoff values is identified
    
        for i=1:10000   %mLS values are calculated for randomly sampled beta and cutoff values below
            cutoff=normrnd(cutoff_mean,cutoff_std);
            beta=normrnd(beta_mean,beta_std);
            constant=y(1,1)/cutoff^beta;
            fit_y=constant*x1.^beta; 
            midx=10^((log10(max(Area))+(log10(cutoff)))/2);
            ac=Refmidy/(11111*Refmidx^beta);
            mLS_array(i,1)=log10((midy/(ac*midx^(beta))));
        end

    mLS_array=mLS_array(all(~isinf(mLS_array),2),:); % "Inf" cells are removed from the array  
    error=std(mLS_array(:));    %Uncertainty of mLS calcultated as a first standard deviation of mLS values
else
    disp('Uncertainty in mLS will not be calculated because the variable "cutoff_error" and "beta_error" is missing')
    error='?'
end

% A graph showing the frequency-area distribution of the given landslides 
% and the corresponding power-law fit are plotted.
loglog(x1,fit_y_stored,'-','LineWidth',2,'Color','r');hold on
loglog(x1,FD,'ok','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','k')
axis([1 1.E+7 1.E-6 1000])
set(get(gca,'Xlabel'),'string','Landslide Area (m^2)','FontSize',12, 'FontUnits','points','FontWeight','normal')
set(get(gca,'Ylabel'),'string','Frequency Density (m^-^2)','FontSize',12, 'FontUnits','points','FontWeight','normal')
str={['\beta = ',num2str(beta_stored)];['mLS = ',num2str(mLS_stored),(char(177)),num2str(error)]};
text(x1(1,1),(min(FD(FD>0)*10)),str,'FontSize',12)  
end
