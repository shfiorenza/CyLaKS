function image2D = imageGaussianKymograph(dataMatrix,siteLength,pixelLength,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot)
%IMAGEGAUSSIANKYMOGRAPH generates a 2D image of a kymograph made from a
%matrix of point positions
%
%SYNOPSIS image2D = imageGaussianKymograph(numPixels,gaussSigma,gaussCenter,gaussAmp,bkgLevel,noiseStd,doPlot)
%
%INPUT  dataMatrix : Motor occupancy matrix, row is occupancy, time increases down
%       siteLength : length (in nm) of single site in occupancy matrix
%       pixelLength: length (in nm) of image pixels
%       gaussSigma : Gaussian standard deviation (pixels, but does not have to be integer).
%       gaussAmp   : Vector of Gaussian amplitudes above background (a.u.).
%       bkgLevel   : Background intensity (a.u.).
%       noiseStd   : Noise standard deviation (a.u.).
%       doPlot     : 1 to plot image, 0 otherwise.
%                    Code both makes a line plot and displays the image
%                    using imshow.
%                    Optional. Default: 0.
%
%OUTPUT image      : 2D image of kymograph
%
%Meredith Betterton, November 2015
%based on imageGaussians1D by Khuloud Jaqaman, October 2011

%% Input

%check number of input arguments
if nargin < 7
    error('imageGaussianAster2D: Missing input arguments')
end
if nargin < 8 || isempty(doPlot)
    doPlot = 0;
end

%% Output

numPixelsX = ceil(size(dataMatrix,2)*siteLength/pixelLength);
numPixelsY = size(dataMatrix,1);

%generate image of background level and noise
image2D = bkgLevel*ones(numPixelsY,numPixelsX) + randn(numPixelsY,numPixelsX)*noiseStd;


%% Generate image

%define sqrt(2pi) to avoid calculating it over and over again
sqrt2Pi = sqrt(2*pi);

%get coordinates of left and right pixel edges
coordLeftEdge = (0.5:numPixelsX-0.5);
coordRightEdge = (1.5:numPixelsX+0.5);

for jTime = 1:numPixelsY     %loop over time points
    image1D = zeros(1,numPixelsX);
    
    occupancyVector = dataMatrix(jTime,:);  %extract the occupancy data
    numGauss = nnz(occupancyVector); 
    gaussCenter = find(occupancyVector) * siteLength/pixelLength;
    gaussAmpVec = ones(size(gaussCenter)) * gaussAmp;
    
    %determine the contribution of each Gaussian to the image
    for iGauss = 1 : numGauss
        
        %calculate the CDF at the left and right pixel edges
        gaussCDFLeftEdge = normcdf(coordLeftEdge,gaussCenter(iGauss),gaussSigma);
        gaussCDFRightEdge = normcdf(coordRightEdge,gaussCenter(iGauss),gaussSigma);
        
        %thus calculate the intensity in each pixel from this Gaussian
        pixelIntensity = (gaussCDFRightEdge - gaussCDFLeftEdge) ...
            * gaussAmpVec(iGauss) * gaussSigma * sqrt2Pi;
                
        image1D = image1D + pixelIntensity;
    end
    
    %add to appropriate row in the kymograph
    image2D(jTime,:) = image2D(jTime,:) + image1D;
end


%% Plot image

if doPlot
    figure, imagesc(image2D); 
    colormap gray; 
    %axis equal
end

end

