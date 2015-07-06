using PureSeq
using Base.Test


## BamReader

# test forward reads
reader = BamReader("data/small.bam", :forward, ReferenceContigs_hg38)
@test position(reader) == 10544
advance!(reader)
@test position(reader) == 11247
close(reader)

# test reverse reads
reader = BamReader("data/small.bam", :reverse, ReferenceContigs_hg38)
@test position(reader) == 10056
advance!(reader)
@test position(reader) == 10102
close(reader)

# test reading through a whole file
reader = BamReader("data/small.bam", :any, ReferenceContigs_hg38)
lastPos = -1
while !eof(reader)
	lastPos = position(reader)
	advance!(reader)
end
close(reader)

# test reading through a whole file with an iterator
reader = BamReader("data/small.bam", :reverse, ReferenceContigs_hg38)
for pos in eachposition(reader)
end
close(reader)


## FeatureMap

# test forward reads
reader = BamReader("data/small.bam", :forward, ReferenceContigs_hg38)
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
fm1 = FeatureMap(BamReader("data/small.bam", :forward, ReferenceContigs_hg38), ones(21)*2)
fm2 = FeatureMap(BamReader("data/small.bam", :forward, ReferenceContigs_hg38), ones(21))
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
fm1 = FeatureMap(BamReader("data/small.bam", :forward, ReferenceContigs_hg38), ones(21)*2)
fm2 = FeatureMap(BamReader("data/small.bam", :forward, ReferenceContigs_hg38), ones(21))
count = 0
for block in denseblocks([fm1, fm2], 1000, loop=true)
	count += 1
	if count > 12000 break end
end
@test count == 12001
close(fm1)
close(fm2)


## SamWriter

f = open("data/tmp.sam", "w")
sw = SamWriter(f, ReferenceContigs_hg38)
writeRead(sw, 10, 16)
close(f)

## write_binned and BinnedReader

# test forward reads
write_binned("data/small.bam", 1000, :reverse, skipDup=false)
reader = BinnedReader("data/small.bam.rbin1000")
@test position(reader) == 11 # where data first enters the window
@test value(reader) == 3.0

while !eof(reader)
	if position(reader) == 12
		@test value(reader) == 5.0
	end
	advance!(reader)
end

# now skip duplicate reads
write_binned("data/small.bam", 1000, :reverse, skipDup=true)
reader = BinnedReader("data/small.bam.rbin1000")
@test position(reader) == 11 # where data first enters the window
@test value(reader) == 3.0
while !eof(reader)
	if position(reader) == 12
		@test value(reader) == 4.0
	end
	advance!(reader)
end
close(reader)


## ContextMap

# test forward reads
reader = BamReader("data/small.bam", :forward, ReferenceContigs_hg38)
cm = ContextMap(reader, 1000, 1000)
@test position(cm) == 9544 # where data first enters the window
@test value(cm) == 1.0
while !eof(cm)
	if position(cm) == 10544
		@test value(cm) == 3.0
	end
	advance!(cm)
end
close(cm)

# test bins
reader = BinnedReader("data/small.bam.fbin1000")
cm = ContextMap(reader, 2, 2)
@test position(cm) == 9 # where data first enters the window
@test value(cm) == 1.0
while !eof(cm)
	if position(cm) == 12
		@test value(cm) == 13.0
	end
	advance!(cm)
end
close(cm)
