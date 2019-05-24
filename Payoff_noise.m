function score=Payoff_noise(stra_player1,stra_player2,b,thred)  
 
 switch (2*stra_player1+stra_player2)
      case 0,score=1;
      case 1,score=-0.15;
      case 2,score=b;
      case 3,score=0.04;
 end
 score=score+unifrnd(0,thred); % Add some noise in the score
%  score=score+normrnd(0,thred); % Add some noise in the score

end