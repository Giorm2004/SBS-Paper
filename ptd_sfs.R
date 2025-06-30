library("pracma")
library("ptdalgorithms")
#inputs: k, k1, k2, m12, m21, IPV sorted in decreasing order of lines in pop 1, with all non-singleton states excluded
sfs_graph <- function(k, k1, k2, m12, m21, IPV) {

#getting fictional states, need to make sure every state actually has proper number of lines
  graph <- create_graph(2 * k)

  start <- vertex_at(graph, 1)

  absorb_state <- rep(0, 2*k)
absorb_state[2*k] <- 1
absorb <- find_or_create_vertex(graph, absorb_state)

initial <- rep(0, 2*k)
initial[1] <- k
find_or_create_vertex(graph, initial)

index <- 3

while (index <= vertices_length(graph)) {
	vertex <- vertex_at(graph, index)
	state <- vertex$state
	#print(index)
	#print(state)
	#check if possible initial state
	if (state[1] + state[k+1] == k) {
		#print("possible initial state")
		#print(IPV[state[k+1] + 1])
	add_edge(start, find_or_create_vertex(graph, state), IPV[state[k+1] + 1])
	#print("initial state added")
	}
	for (i in 1:(2*k)){
	#set rates for pop 1s
	if ((i %in% 1:(k-1)) & (state[i] > 0)){ # nolint: vector_logic_linter.
		#coalescence of two branches subtending same number of lineages
		if (state[i] > 1){
			if (i * 2 == k){
				add_edge(vertex, absorb, 1/k1)
			}
			else{
		#print("equal line coal")
		child <- sapply(state, function(i) i)
		child[i] <- state[i] - 2
		child[2*i] <- state[2*i] + 1
		#print(child)
		add_edge(vertex, find_or_create_vertex(graph, child), nchoosek(state[i], 2)/k1)
			}}
		#migration
		#print("migration")
		child <- sapply(state, function(i) i)
		child[i] <- state[i] - 1
		child[i+k] <- state[i+k] + 1
		#print(child)
		add_edge(vertex, find_or_create_vertex(graph, child), state[i]*m12)
		#coalescence of branches subtending different number of lineages
		for (j in (i+1):(k-1)){
		if ((state[j] >= 1) & (i + j <= k)){ # nolint: vector_logic_linter.
			
		if (i + j == k) {
		#absorbing
		#print("absorbing")

		add_edge(vertex, absorb, 1/k1)
		}
		else{
			#print("unequal line coal")

		child <- sapply(state, function(i) i)
		child[i] <- state[i] - 1
		child[j] <- state[j] - 1
		child[i+j] <- state[i+j] + 1
		#print(child)
		add_edge(vertex, find_or_create_vertex(graph, child), state[i]*state[j]/k1)
		}
		


		}
		}
	
	}

	if ((i %in% (k+1):(2*k-1)) & (state[i] > 0)){ # nolint: vector_logic_linter.
		#print('down here now')
	#coalescence of two branches subtending same number of lineages
		if (state[i] > 1){
			#print("equal line coal")
			if (2*i == 3 * k){
				add_edge(vertex, absorb, 1/k2)
			}
			else{
	child <- sapply(state, function(i) i)
	child[i] <- state[i] - 2
	child[2 * i - k] <- state[2 * i - k] + 1
	#print(child)
	add_edge(vertex, find_or_create_vertex(graph, child), nchoosek(state[i], 2)/k2)
			}}
	#migration
	#print('migration')
	child <- sapply(state, function(i) i)
	child[i] <- state[i] - 1
	child[i-k] <- state[i-k] + 1
	#print(child)
	add_edge(vertex, find_or_create_vertex(graph, child), state[i]*m21)
#coalescence of branches subtending different number of lineages
for (j in (i+1):(2*k-1)){
if ((state[j] >= 1) & (i + j <= 3*k)){ # nolint: vector_logic_linter.
if (i + j == 3*k) {
#absorbing
#print('absorb')
add_edge(vertex, absorb, 1/k2)
}
else{
#print('unequal line coal')
#print(c(i, j))
child <- sapply(state, function(i) i)
child[i] <- state[i] - 1
child[j] <- state[j] - 1
child[i+j - k] <- state[i+j - k] + 1
#print(child)
add_edge(vertex, find_or_create_vertex(graph, child), state[i]*state[j]/k2)
}



}
}


}

	}
#print(edges(vertex))
index <- index + 1
}

#print("finished")
return (graph)
}
binom_prob <- function(k,n,p){
	return (nchoosek(n, k) * (p**k) * ((1-p)**(n-k)))
}
binom_vec <- function(n,p){
	return (sapply(0:n, binom_prob, n = n, p = p))

}

ipv <- binom_vec(10, 0.5)
start <- proc.time()[3]
gr <- sfs_graph(20, 0.25, 0.25, 0.054,  0.05, ipv)
end <- proc.time()[3]

sbs_graph <- function(k, p, ph, rho){
  sfs_graph(k, ph, ph, rho / 2, rho / 2, binom_vec(k, p))
}

overdom_graph <- function(k, p, rho){
	sfs_graph(k, p, 1-p, rho * (1-p), rho*p, binom_vec(k, p))
}

nsfs <- function(k){
    sapply(1:(k-1), function(i) 1/i)
}
marg_moments <- function(i, k, graph, ord) {
  rewards <- states(graph)[, i] + states(graph)[, i + k]
  moments(graph, ord, rewards)
}
mm_spec <- function(k, graph, ord) {
  sapply(1:(k-1), marg_moments, k = k, graph = graph, ord = ord)
}
print(sapply(1:9, marg_moments, k = 10, graph = gr, ord = 1 ))
fold <- function(i, usfs){
  l <- length(usfs)
  if (2 * i == l){
    return(usfs[i])
  
  usfs[i] + usfs[l-i]
}
}
foldsfs <- function(usfs){
  l <- length(usfs) %/% 2
  sapply(1:l, fold, usfs = usfs)
}
print(end - start)