q.update <-
function(z){
  rbeta(1, 1+sum(z==1), 1 + sum(z==2))
}
