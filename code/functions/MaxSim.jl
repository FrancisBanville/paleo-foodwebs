"""
MaxSim(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns the average of all species’ largest similarity index 
"""
function MaxSim(N::UnipartiteNetwork)
      # additive Jaccard similarity between all species pairs
      AJS_N = AJS(N) 
      # find the maximum similarity for every species
      length_AJS_N = length(AJS_N)
      max_AJS = []
            for i in 1:richness(N) 
                  spi = []
                  for j in 1:length_AJS_N
                        if in(AJS_N[j][1])(species(N)[i])
                              push!(spi, AJS_N[j][2])
                        end
                  end
                  if length(spi) > 0
                        push!(max_AJS, maximum(spi)[1])
                  end
            end
      # return average similarity index 
      return mean(max_AJS) 
end
