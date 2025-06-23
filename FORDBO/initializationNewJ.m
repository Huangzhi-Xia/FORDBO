%% Good nodes set method
function Positions=initializationNew(SearchAgents_no,dim,ub,lb)


PositionsG = Goodnode(SearchAgents_no,dim);
Positions = PositionsG.*(ub-lb)+lb;

end




