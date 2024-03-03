%% analytic estimate the threshold based on Theorem 1

function thre = tail_est ( alpha, M )


 for b = 2 : 0.01 : 4;

    Pr=0;

    for B = 2:M
    
        temp1 = exp(-0.5*(b^2));
        temp2 = b*sqrt((2*B-1)/B/(B-1));
        temp2 = fv(temp2);
        temp3 = (b^2)*(2*B-1)/2/B/(B-1)/sqrt(2* pi);
        Pr = Pr + temp1 * temp2 * temp3;
    
    end

 
   
   if Pr <= alpha 
     
      thre = b;
      
      break
   
   end
     
   
 end


end
