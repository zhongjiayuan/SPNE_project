function p = SNN(deltpcc,pre_pcc,num)
if abs(pre_pcc)==1
    pre_pcc=0.99999999;
end
 z=deltpcc/((1-pre_pcc^2)/(num-1));
 p=1-normcdf(abs(z));
end