using PureSeq
using Base.Test

## BamReader

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


## FeatureMap

# test forward reads
reader = BamReader("data/small.bam", false, ReferenceContigs_hg38)
fm = FeatureMap(reader, ones(2001))
@test position(fm) == 9544 # where data first enters the window
@test value(fm) == 1.0
while !eof(fm)
	if position(fm) == 10544
		@test value(fm) == 3.0
	end
	advance!(fm)
end
close(fm)

## DenseBlockIterator

# test forward reads
fm1 = FeatureMap(BamReader("data/small.bam", false, ReferenceContigs_hg38), ones(21)*2)
fm2 = FeatureMap(BamReader("data/small.bam", false, ReferenceContigs_hg38), ones(21))
count = 0
for block in denseblocks([fm1, fm2], 100)
	if count == 105
		@test block[33,1] == 0.0
		@test block[33,2] == 0.0
		@test block[34,1] == 2.0
		@test block[34,2] == 1.0
		@test block[54,1] == 2.0
		@test block[54,2] == 1.0
		@test block[55,1] == 0.0
		@test block[55,2] == 0.0
	end
	count += 1
	if count > 120 break end
end
close(fm1)
close(fm2)

# test infinite looping of blocks
fm1 = FeatureMap(BamReader("data/small.bam", false, ReferenceContigs_hg38), ones(21)*2)
fm2 = FeatureMap(BamReader("data/small.bam", false, ReferenceContigs_hg38), ones(21))
count = 0
for block in denseblocks([fm1, fm2], 1000, loop=true)
	count += 1
	if count > 12000 break end
end
@test count == 12001
close(fm1)
close(fm2)