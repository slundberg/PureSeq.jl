using PureSeq
using Base.Test

# test forward reads
reader = BamReader("data/small.bam", false, ReferenceContigs_hg38)
@test position(reader) == 10544
advance!(reader)
@test position(reader) == 11247
close(reader)

# test reverse reads
reader = BamReader("data/small.bam", true, ReferenceContigs_hg38)
@test position(reader) == 10056
advance!(reader)
@test position(reader) == 10102
close(reader)

# test reading through a whole file
reader = BamReader("data/small.bam", true, ReferenceContigs_hg38)
while !eof(reader)
	advance!(reader)
end
close(reader)

# test reading through a whole file with an iterator
reader = BamReader("data/small.bam", true, ReferenceContigs_hg38)
for pos in eachposition(reader)
end
close(reader)