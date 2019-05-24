function score=Payoff(stra_player1,stra_player2,b)  
 
 switch (2*stra_player1+stra_player2)
      case 0,score=1;
      case 1,score=-0.15;   %-0.15;
      case 2,score=b;
      case 3,score=0.04;   %0.04;
 end
 %score=score+normrnd(0,0.04); % Add some noise in the score
end