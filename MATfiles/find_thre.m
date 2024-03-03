%% analytic computation
function thre = find_thre(M, ARL)


 for b=1.5:0.01:5

    temp1=exp(-0.5*(b^2));
    temp2=1/sqrt(2* pi);
    temp3=b*sqrt(2*(2*M-1)/M/(M-1));
    temp4=fv(temp3);
    temp5=(b^2)*(2*M-1)/M/(M-1);
    temp6=temp4*temp5;
 
    Pr=temp6*temp1*temp2;
    
    
    if ( 1/Pr >= ARL )
      
      thre = b;
      break;
      
    end


 end

end
