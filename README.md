# vcf-filter
Java tool to filter an input VCF file a la carte, based on HTSJDK and Apache commons CLI.

# dependencies
|Software|Source|
|--------|------|
|Java 1.8|<https://www.java.com/es/download>|
|picard.jar , included in lib/|<https://github.com/broadinstitute/picard>|
|commons-cli-1.3.1.jar, included in lib/|<https://commons.apache.org/proper/commons-cli>|

# compilation instructions
javac -cp lib/picard.jar:lib/commons-cli-1.3.1.jar:. Split.java VcfReader.java

# check supported operations
If your system has Perl installed, as most Linux/Mac OS X systems, you can simply call:
$ ./vcf-filter 

Alternatively, just call it directly with java:
$ java -cp lib/picard.jar:lib/commons-cli-1.3.1.jar:. VcfReader

# author
Eduardo Candeal, with help from Carlos P Cantalapiedra and Bruno Contreras-Moreira
