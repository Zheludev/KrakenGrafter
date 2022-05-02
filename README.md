# KrakenGrafter
a simple python3 script for "grafting" on novel sequences to a custom [Kraken2](https://github.com/DerrickWood/kraken2) database

the purpose of this script is to append on the contents of a `.fasta` file to the `nodes.dmp` and `names.dmp` used by [Kraken2](https://github.com/DerrickWood/kraken2) - in effect "grafting" on a new branch to the tree of life inside the `nodes.dmp` and `names.dmp` files. this allows the user to then use the `kraken2-build --add-to-library` functionality during database building thereby making sure their own custom sequences are included

this is an extension to simply using `kraken2-build --add-to-library` because it allows the user to create new nodes that did not previously exist in the taxonomic system - see [kraken2's custom database section](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases).

=================================================

two example usecases: 
- say you've identified a new set of members of a genus and want to count their prevalance in sequencing datasets - they are hitherto unknown to the NCBI, and so wont have their own nodes in the taxonomy - you'd use the ["A)" approach here](https://github.com/Zheludev/KrakenGrafter/edit/main/README.md#1-a-given-fasta-files-sequences-isare-inserted-directly-beneath-taxonomically-speaking-a-specified-root-node)
- say you've indetified a completely new set of sequences - not only does the NCBI doesn't have a node, but maybe they don't even fit the tree of life! - you could append these sequences to "synthetic construct" under a new genus-species pair - you'd use the ["B)" approach here](https://github.com/Zheludev/KrakenGrafter/edit/main/README.md#2-the-user-declares-a-brand-new-parent-node-that-is-inserted-beneath-the-root-node-at-a-user-specified-taxonomic-depth-and-then-the-given-fasta-sequences-isare-inserted-directly-beneath-this-new-parent-node-this-functionality-is-enabled-by-declaring-both--parent_taxon-and--parent_rank)

=================================================

**installation**

download the KrakenGrafter.py file and place it somewhere in your `PATH` - make sure you are using python3 with [argparse](https://docs.python.org/3/library/argparse.html) installed

=================================================

**the script takes in three input files:**
* `nodes.dmp` to modify
* `names.dmp` to modify
* `.fasta` of sequences to graft

NOTE: the seqIDs need to contain no whitespace and no pipe characters ("|")

the entire seqID will be grafted into the `nodes.dmp` and `names.dmp` files so maybe keep it concise

=================================================

**the script then generates three output files:**
* `new_nodes.dmp` (or specified name) of the modified `nodes.dmp`
* `new_names.dmp` (or specified name) of the modified `names.dmp`
* `K2.fasta` (or specified name) of the modified `.fasta` file

=================================================

the script can be used in two different ways - in both cases, brand new taxon IDs are generated that have not been seen in the `nodes.dmp` and `names.dmp` files:
#### A) a given `.fasta` file's sequence(s) is/are inserted directly beneath (taxonomically speaking) a specified "root" node

  example:

  `python KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta input.fasta -root 32630`

=================================================

#### B) the user declares a brand new "parent" node that is inserted beneath the "root" node (at a user-specified taxonomic depth) and then the given `.fasta` sequence(s) is/are inserted directly beneath this new "parent" node this functionality is enabled by declaring both `-parent_taxon` and `-parent_rank`

  example:

  `python KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta input.fasta -root 32630 -parent_taxon new_genus -parent_rank genus`

=================================================

NOTE: `i_nodes`, `i_names`, `root`, `o_nodes`, `o_names`, `o_fasta`, and `debug` are all optional variables to declare.

the default "root" node is 32630 - which is the NCBI's "synthetic construct" node - a safe habour for new sequences.

To find a more specific node, one can either go on the [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) or use `grep 'taxon_string' names.dmp` to try and work out what taxon ID to use for "root"

=================================================

=================================================

## tutorial

regardless of how you use `KrakenGrafter.py`, you need to start the same way by downloading the NCBI taxonomy:

1. `kraken2-build --download-taxonomy --db $DB_name`

    then, optionally, you can download any number of the [prexisting kraken2 databases](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases) e.g.:

    `kraken2-build --download-library viral --db $DB_name`

    `kraken2-build --download-library plant --db $DB_name`
    
    if you are having issues with the FTP path when running `kraken2-build --download-library` then [this is a likely fix](https://github.com/DerrickWood/kraken2/issues/508) (in `rsync_from_ncbi.pl` replace `if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##)) {` with `if (! ($full_path =~ s#^https://${qm_server}${qm_server_path}/##)) {`)

2. now that you have the databases of interest (please see [this bioRxiv](https://bit.ly/3EWkYJf) for an assessment of how to choose your database and classification settings), navigate to:

    `cd $DB_name/taxonomy`

    here you will find the `nodes.dmp` and `names.dmp` files needed for `KrakenGrafter.py`

    I suggest you copy and re-name these original `nodes.dmp` and `names.dmp` files and keep them in `$DB_name/taxonomy` - that way you can just revert their names and use them as intended if something goes wrong.

3. you also need to identify the taxon ID (taxID) under which the grafting will be performed, e.g.:

    `grep 'Enterobacterales' names.dmp` _(the first column will contain the taxID - an integer - be sure to manually inspect)_

4. lastly in the prep, you need to have ready a `.fasta` file. The file can contain any number of sequences, the only thing to look out for are the seqIDs (the ">" names of the sequences) - **your `.fasta` file seqIDs need to have no whitespaces (space / tab) and no pipe characters ("|")**

    **your `.fasta` file must contain sequences that will all be insterted under the same taxonomic node** - e.g. if you want to graft on both say a new bacterium and a new protist, you would run `KrakenGrafter.py` twice, once for the bacterium and once for the protist (any order) - each with their own `.fasta` file - _but the input `nodes.dmp` and `names.dmp` files of the second 'round' of `KrakenGrafter.py` will be the output from the first_

5. I'd now take your copies of the `nodes.dmp` and `names.dmp` files and your `.fasta` file(s) and put them in the same directory (this isn't necessary but I find it easier than writing out long paths)

### A) grafting on new sets of sequences to pre-existing nodes

   let's say I want to graft on `species1.fasta` under `genus_taxID` and `genus2.fasta` under `family_taxID` - this would be done in two steps:
   
 1) `python KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta species1.fasta -root genus_taxID -o_nodes nodes.dmp -o_names names.dmp -o_fasta K2_species1.fasta`
   
    this overwrites the input `nodes.dmp` and `names.dmp` files
   
 2) `python KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta genus2.fasta -root family_taxID -o_nodes nodes.dmp -o_names names.dmp -o_fasta K2_genus2.fasta`
   
    note how in this set-up, the output `nodes.dmp` and `names.dmp` files of 1) are used in 2) - that way all the sequences are added to same `nodes.dmp` and `names.dmp` files
   
 3) **now make sure to move the edited `nodes.dmp` and `names.dmp` files (with those exact names) back to `$DB_name/taxonomy`**
   
 4) now we need to use the [in-built features of Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases) to actually add these sequences to the database (we can do this in one step by catting the `KrakenGrafter.py` output `.fasta` files) and build the database:
   
    `cat K2_species1.fasta K2_genus2.fasta > K2_seqs_to_graft.fasta`
   
    `kraken2-build --add-to-library K2_seqs_to_graft.fasta --db $DB_name`
   
    `kraken2-build --build --db $DB_name --threads 32`

### B) grafting on new sets of sequences to new nodes

   let's say I want to graft on `species3.fasta` to a new node under a given `family_taxID` in a new genus (named `new_genus`) - in addition to the setup above, you will also need to decide on a name for the new parent taxon (e.g. `new_genus`) as well as a taxonomic rank for this taxon (in this case `genus`). The 'root' under which this parent taxon is inserted can be any taxonomic depth above the parent taxon (e.g. you could instert a new bacterial species in simply the bacterial kingtom taxID as the 'root' - 2 - but then create a brand new genus - `new_genus`). Because the new parent taxon (in this case `new_genus`) doesn't yet exist in the taxomony, a new taxID will be automatically generated for it as well as for the actual sequences being added below it:
 
   1) `python KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta species3.fasta -root family_taxID -parent_taxon new_genus -parent_rank genus -o_nodes nodes.dmp -o_names names.dmp -o_fasta K2_species3.fasta`

   2) **now make sure to move the edited `nodes.dmp` and `names.dmp` files (with those exact names) back to `$DB_name/taxonomy`**
   
   3) now we need to use the [in-built features of Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases) to actually add these sequences to the database and build the database:
      
    `kraken2-build --add-to-library K2_species3.fasta --db $DB_name`
   
    `kraken2-build --build --db $DB_name --threads 32`

=================================================

written by INZ - 04/21/22, Stanford Unversity, provided with no acceptance of liability or promise of functionality, version 0.1.0

=================================================
