# MarkDuplicates Rust
- A Rust porting version of Picard MarkDuplicates
- Made for testing there is any performance benefit when writing it with Rust.

### 🔭 This is under development.
- Currently implemented up to writing output bam.
- Metrics writing code should be written.


### 🐃 Interim results
I tested performance with NA12878.mapped.bam, I was a bit disappointed with the results.  
<br>
MarkDuplicates actually consumed most of time(About 70%) to write bam.  
With 1 thread, rust-htslib and htsjdk have almost the same speed for writing bam.  
<br>
Finding duplicates of rust implementation was about 1.6x faster than java(original Markduplicates), but this took up only a small fraction of the total running time.  
<br>
Using 4 or 8 threads for writing and reading bam, rust implementation was much faster than java (11min, 8min vs 29min) but if htsjdk support multiple thread for writing and reading bam, I think the performance benefit may be little.



