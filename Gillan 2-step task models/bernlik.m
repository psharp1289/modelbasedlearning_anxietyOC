function prob=bernlik(N,P,hits)
prob=P^hits*(1-P)^(N-hits);
end