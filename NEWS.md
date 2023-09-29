# pfGO 1.2

* Formatting fixes to all significant-genes-in-significant-terms output-files for easier downstream parsing.
* Functions pulling data from PlasmoDB.org updated to 'current' release
    * no longer hard-coded to specific version release--should work to access all "Current_Release" links so long as there are no updates to PlasmoDB url structure or file-compression (works with non-compressed files and gzip-compressed files)
    * get.annot function now automatically pulls all 'gene' records from .gff files ("protein_coding_gene","ncRNA_gene","pseudogene")
* Included example-datasets updated from latest PlasmoDB release (version 66)


# pfGO 1.0

* Added a `NEWS.md` file to track changes to the package.
