export BinningMap, close, value, position, eof, advance!

type BinningMap
    reader::BamReader
    binSize::Int64
    position::Int64
    value::Float64
end

function BinningMap(reader::BamReader, binSize)
    fm = BinningMap(reader, binSize, 0, 0.0)
    
    advance!(fm)
    fm
end
close(fm::BinningMap) = close(fm.reader)
value(fm::BinningMap) = fm.value
position(fm::BinningMap) = fm.position
eof(fm::BinningMap) = fm.position <= 0

function advance!(fm::BinningMap)
    fm.position = floor((fm.reader.position-1)/fm.binSize) + 1
    binEnd = fm.position*fm.binSize
    
    # Fill in the bin
    fm.value = 0.0
    while fm.reader.position != -1 && fm.reader.position <= binEnd
        fm.value += 1
        PureSeq.advance!(fm.reader)
    end
end