function image2D = imageGaussianOverlap(dataVector,siteLength,pixelLength,...
    pixelPad,gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot)
%IMAGEGAUSSIANOVERLAP generates a 2D image of a motor distribution made
%from an occupancy vector
%
%SYNOPSIS image2D = imageGaussianOverlap

%INPUT  dataVector : Motor occupancy vector, row vector
%       siteLength : length (in nm) of single site in occupancy matrix
%       pixelLength: length (in nm) of image pixels
%       pixelPad   : pixel padding at ends of overlap and in y
%       gaussSigma : Gaussian standard deviation (pixels, but does not have to be integer).
%       gaussAmp   : Vector of Gaussian amplitudes above background (a.u.).
%       bkgLevel   : Background intensity (a.u.).
%       noiseStd   : Noise standard deviation (a.u.).
%       doPlot     : 1 to plot image, 0 otherwise.

%OUTPUT image      : 2D image of overlap
%
%Meredith Betterton, November 2015
%based on imageGaussians1D by Khuloud Jaqaman, October 2011

%% Input

%check number of input arguments
if nargin < 8
    error('imageGaussianAster2D: Missing input arguments')
end
if nargin < 9 || isempty(doPlot)
    doPlot = 0;
end

%% Output

numPixelsX = ceil(size(dataVector,2)*siteLength/pixelLength)+2*pixelPad;
numPixelsY = 2*pixelPad;

%generate image of background level
image2D = bkgLevel*ones(numPixelsY,numPixelsX);

%% Generate image

%get coordinates of left and right pixel edges
coordLeftEdge = (0.5:numPixelsX-0.5);
coordRightEdge = (1.5:numPixelsX+0.5);
coordTopEdge = (0.5:numPixelsY-0.5)';
coordBottomEdge = (1.5:numPixelsY+0.5)';


occupancyVector = dataVector;  %extract the occupancy data
numGauss = nnz(occupancyVector);
gaussCenter = find(occupancyVector) * siteLength/pixelLength + pixelPad;
gaussAmpVec = occupancyVector(find(occupancyVector)) * gaussAmp;

%determine the contribution of each Gaussian to the image
for iGauss = 1 : numGauss
    
    %calculate the CDF at the left and right pixel edges
    gaussCDFLeftEdge = normcdf(coordLeftEdge,gaussCenter(iGauss),gaussSigma);
    gaussCDFRightEdge = normcdf(coordRightEdge,gaussCenter(iGauss),gaussSigma);

    %calculate the CDF at the top and bottom pixel edges
    gaussCDFTopEdge = normcdf(coordTopEdge,pixelPad,gaussSigma);
    gaussCDFBottomEdge = normcdf(coordBottomEdge,pixelPad,gaussSigma);

    
    %thus calculate the intensity in each pixel from this Gaussian
    pixelIntensity = (gaussCDFBottomEdge - gaussCDFTopEdge) *...
        (gaussCDFRightEdge - gaussCDFLeftEdge) *...
        gaussAmpVec(iGauss) /(2*pi);
    
    image2D = image2D + pixelIntensity;
end
%add noise
image2D = image2D + (image2D>(2*bkgLevel)) .* randn(numPixelsY,numPixelsX)*noiseStd;


%% Plot image

if doPlot
    figure, imagesc(image2D); 
    colormap gray; 
    %axis equal
end

end

