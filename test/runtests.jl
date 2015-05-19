using PureSeq
using Base.Test

# test forward reads
reader = BamReader("data/small.bam", false, ReferenceContigs_hg38)
@test reader.position == 10544
advance!(reader)
@test reader.position == 11247
close(reader)

# test reverse reads
reader = BamReader("data/small.bam", true, ReferenceContigs_hg38)
@test reader.position == 10056
advance!(reader)
@test reader.position == 10102
close(reader)

# test reading through a whole file
reader = BamReader("data/small.bam", true, ReferenceContigs_hg38)
while reader.position != -1
	advance!(reader)
end
close(reader)