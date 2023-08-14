"""
MaxSim(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns the average of all species’ largest similarity index 
"""
function MaxSim(N::UnipartiteNetwork)
      # additive Jaccard similarity between all species pairs
      AJS_N = AJS(N) 
      
      # convert to vectors to facilitate computations
      n = length(AJS_N)
      sim = [AJS_N[i][2] for i in 1:n] 
      sp1 = [collect(AJS_N[i][1])[1] for i in 1:n]
      sp2 = [collect(AJS_N[i][1])[2] for i in 1:n]
    
      # find the maximum similarity for every species
      sp = species(N)
      MxSim_sp = zeros(Float64, length(sp))
      for (i, sp) in enumerate(species(N))
            MxSim_sp[i] = maximum(vcat(sim[sp1 .== sp], sim[sp2 .== sp]), init=0)
      end
      
      # return average similarity index (of all species except those without similarity values)
      return mean(MxSim_sp[Not(MxSim_sp .== 0)])
end



