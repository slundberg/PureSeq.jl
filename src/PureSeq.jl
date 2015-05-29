module PureSeq

import Base: eof

include("ReferenceContigs.jl")
include("BamReader.jl")
include("FeatureMap.jl")
include("DenseBlockIterator.jl")
include("SamWriter.jl")
include("BinningMap.jl")

end