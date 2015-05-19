export BamReader, close, advance!

type BamReader
    bamStream
    useReverseReads::Bool
    done::Bool
    position::Int64
    contigs::ReferenceContigs
end

function BamReader(bamFileName::ASCIIString, useReverseReads, contigs)
    f = GZip.open(bamFileName)
    
    # make sure this is a BAM file
    code = read(f, Uint8, 4)
    @assert code == b"BAM\1"
    
    # get through the header data 
    l_text = read(f, Int32)
    skip(f, l_text)
    
    # make sure the contigs match our reference
    n_ref = read(f, Int32)
    @assert n_ref == contigs.count
    for j in 1:n_ref
        l_name = read(f, Int32)
        refName = convert(ASCIIString, read(f, Uint8, l_name)[1:end-1]) # ignore the null terminator
        l_ref = read(f, Int32)
        @assert l_ref == contigs.sizes[j]
        @assert refName == contigs.names[j]
    end
    
    r = BamReader(f, useReverseReads, false, 1, contigs)
    advance!(r)
    r
end
function close(reader::BamReader)
    GZip.close(reader.bamStream)
end

function advance!(r::BamReader)
    f = r.bamStream
    while !r.done
        if peek(f) == -1 # eof does not work on the BAM files either in C++ or here (BAM vs. gzip issue?)
            r.done = true
            r.position = -1
            return
        end
        
        buf = Array(Int32, 5) # [block_size, refID, pos, bin_mq_nl, flag_nc]
        gzread(f, pointer(buf), 20)
        block_size = buf[1]
        refID = buf[2] + 1 # the reference contig this read maps to
        
        # get the read position
        if refID != 0
            r.position = buf[3] + r.contigs.offsets[refID] + 1; # convert to 1 based indexing
        end
        
        forward = (buf[5] & 1048576) == 0 # see if we are reverse complemented
        skip(f, block_size-16) # skip the rest of the entry
        
        # break if we found a read in the right direction
        if refID != 0 && ((forward && !r.useReverseReads) || (!forward && r.useReverseReads))
            return
        end
    end
end