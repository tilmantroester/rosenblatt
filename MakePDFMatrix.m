function PDF = MakePDFMatrix

x = linspace(-4,6,1000);
D = [.05:.05:.45];
PDF = zeros(11,1000);
m  = x > -1/sqrt(2);
PDF(1,m) = 1./sqrt(pi + sqrt(2)*pi*x(m)) .* exp(-.5 - x(m)/sqrt(2) ); 
PDF(end,:) = normpdf(x);

parfor i=1:9
  PDF(i+1,:) = RosenblattCDF(x,D(i));
end

save PDF_data


end 
