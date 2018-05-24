function [sd,sm] = integratedsimilarity(FS,FSP,SS,SSP,kd,km)         
sm = FS.*FSP+km.*(-(FSP-1));            
sd = SS.*SSP+kd.*(-(SSP-1));            
end