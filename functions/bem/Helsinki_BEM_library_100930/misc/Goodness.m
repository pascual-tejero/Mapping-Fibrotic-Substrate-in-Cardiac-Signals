function result=Goodness(phi_meas,phi_calc)
%function result=Goodness(phi_meas,phi_calc)
%Goodness of Fit; 1-RE^2
temp=phi_meas-phi_calc;
temp2=temp.^2;
phi_meas2=phi_meas.^2;
rel=sum(temp2)/sum(phi_meas2);
result=1-rel;