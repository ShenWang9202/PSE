function mixtureDistributions(pd1,pd2)
% Derive Cupid distribution objects from MATLAB ones:
ckern1 = dMATLABc(pd1);
ckern2 = dMATLABc(pd2);
% Make a Convolution distribution from the Cupid distribution objects:
convkern = Convolution(ckern1,ckern2);
%convkern = Difference(ckern1,ckern2);

% Plot PDF and CDF
convkern.PlotDens;  
% Compute various properties of the convolution distribution:
meanValue = convkern.Mean
varianceValue = convkern.Variance
interval = zeros(1,2)
interval(1) = convkern.InverseCDF(0.025)
interval(2) = convkern.InverseCDF(0.975)

end
 