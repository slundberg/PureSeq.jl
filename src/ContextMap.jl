using DataStructures

import Base: eof, close
export ContextMap, close, value, position, eof, advance!

type ContextMap
    reader::BamReader
    contextBefore::Int64
    contextAfter::Int64
    position::Int64
    posQueue::Deque{Int64}
    value::Float64
end

function ContextMap(reader::BamReader, contextBefore, contextAfter)
    cm = ContextMap(reader, contextBefore, contextAfter, 0, Deque{Int64}(), 0.0)
    
    advance!(cm)
    cm
end
close(cm::ContextMap) = close(cm.reader)
value(cm::ContextMap) = cm.value
position(cm::ContextMap) = cm.position
eof(cm::ContextMap) = cm.position < 0

function advance!(cm::ContextMap)
    cm.position += 1
    
    # remove old positions that fell outside the sliding window
    while length(cm.posQueue) > 0 && front(cm.posQueue) < cm.position-cm.contextBefore
        shift!(cm.posQueue)
        cm.value -= 1 # assume all positions have value of one!
    end
    
    # skip ahead to the next spot we have data if our queue is empty
    if length(cm.posQueue) == 0
        cm.position = cm.reader.position-cm.contextBefore
    end
    
    # add new positions in the sliding window
    while cm.reader.position != -1 && cm.reader.position <= cm.position+cm.contextAfter
        cm.value += 1 # assume all positions have value of one!
        push!(cm.posQueue, cm.reader.position) 
        advance!(cm.reader)
    end
end