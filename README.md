# KrakenGrafter
a simple python3 script for "grafting" on novel sequences to a custom [Kraken2](https://github.com/DerrickWood/kraken2) database

the purpose of this script is to append on the contents of a `.fasta` file to the `nodes.dmp` and `names.dmp` used by [Kraken2](https://github.com/DerrickWood/kraken2) - in effect "grafting" on a new branch to the tree of life inside the `nodes.dmp` and `names.dmp` files. this allows the user to then use the `kraken2-build --add-to-library` functionality during database building thereby making sure their own custom sequences are included

this is an extension to simply using `kraken2-build --add-to-library` because it allows the user to create new nodes that did not previously exist in the taxonomic system.

an example usecase is: 
- say you've identified a new set of members of a genus and want to count their prevalance in sequencing datasets - they are hitherto unknown to the NCBI, and so wont have their own nodes in the taxonomy - you'd use the "1)" approach here
- say you've indetified a completely new set of sequences - not only does the NCBI doesn't have a node, but maybe they don't even fit the tree of life! - you could append these sequences to "synthetic construct" under a new genus-species pair - you'd use the "2)" approach here

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
* `new_nodes.dmp` (or specified name) of the modified nodes.dmp
* `new_names.dmp` (or specified name) of the modified names.dmp
* `K2.fasta` (or specified name) of the modified .fasta file

=================================================

the script can be used in two different ways - in both cases, brand new taxon IDs are generated that have not been seen in the `nodes.dmp` and `names.dmp` files:
1) a given `.fasta` file's sequence(s) is/are inserted directly beneath (taxonomically speaking) a specified "root" node

  example:

  `python KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta input.fasta -root 32630`

2) the user declares a brand new "parent" node that is inserted beneath the "root" node (at a user-specified taxonomic depth) and then the given `.fasta` sequence(s) is/are inserted directly beneath this new "parent" node this functionality is enabled by declaring both `-parent_taxon` and `-parent_rank`

  example:

  `python KrakenGrafter.py -i_nodes nodes.dmp -i_names names.dmp -i_fasta input.fasta -root 32630 -parent_taxon new_genus -parent_rank genus`

NOTE: `i_nodes`, `i_names`, `root`, `o_nodes`, `o_names`, `o_fasta`, and `debug` are all optional variables to declare the default "root" node is 32630 - which is the NCBI's "synthetic construct" node - a safe habour for new sequences otherwise, one can either go on the NCBI taxonomy browser or use `grep 'taxon_string' names.dmp` to try and work out what taxon ID to use for "root"

=================================================

written by INZ - 04/21/22, Stanford Unversity, provided with no acceptance of liability or promise of functionality, version 0.1.0

=================================================
