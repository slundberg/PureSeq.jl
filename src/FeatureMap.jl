using DataStructures

import Base: eof, close, position
export FeatureMap, close, value, eof, advance!

type FeatureMap
    reader::BamReader
    feature::Array{Float64,1}
    featureWidth::Int64
    position::Int64
    posQueue::Deque{Int64}
    value::Float64
end

function FeatureMap(reader::BamReader, feature)
    fm = FeatureMap(reader, feature, (length(feature)-1)/2, 0, Deque{Int64}(), 0.0)

    advance!(fm)
    fm
end
close(fm::FeatureMap) = close(fm.reader)
value(fm::FeatureMap) = fm.value
position(fm::FeatureMap) = fm.position
eof(fm::FeatureMap) = fm.position < 0

function advance!(fm::FeatureMap)
    fm.position += 1

    # remove old positions that fell outside the sliding window
    while length(fm.posQueue) > 0 && front(fm.posQueue) < fm.position-fm.featureWidth
        shift!(fm.posQueue)
    end

    # skip ahead to the next spot we have data if our queue is empty
    if length(fm.posQueue) == 0
        fm.position = fm.reader.position-fm.featureWidth
    end

    # add new positions in the sliding window
    while fm.reader.position != -1 && fm.reader.position <= fm.position+fm.featureWidth
        push!(fm.posQueue, fm.reader.position)
        advance!(fm.reader)
    end

    # compute the feature value for this position
    fm.value = 0
    for pos in fm.posQueue
        fm.value += fm.feature[pos - fm.position + fm.featureWidth + 1]
    end
end
