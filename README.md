# Tedna

This is the Tedna repository.

It has first been published in [Bioinformatics in 2014](https://academic.oup.com/bioinformatics/article/30/18/2656/2475636).

It also has been used in [Molecular Ecology in 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14969?af=R).

To use Tedna, clone it and type `make`.

Type `./tedna -h` to get some info on the options.

If you have paired-end reads with insert size 300, you can type:

    ./tedna -1 left.fastq -2 right.fastq -k 61 -i 300 -o output.fasta
    
If you have long reads (such as PacBio), you can type:

    ./tedna -1 reads.fastq -k 61 -m 1000 -t 25 -o output.fasta
    
Do no hesitate to [contact me](mailto:matthias.zytnicki@inra.fr) for any further question.
