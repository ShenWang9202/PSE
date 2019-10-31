function plotMenteCarloEach(id,MonteCarloResult,DeterministicResult,ele)
h=histogram(MonteCarloResult);
hold on
title(string(ele) + id)
line([DeterministicResult DeterministicResult],[0  max(h.Values)],'LineWidth',10)
hold off
end


