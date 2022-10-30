# noncanonical_mechanism
The high-throughput pri-miRNA cleavage assays on 262144 variants

Codes for sequencing data analyses:

0.make_ref.py: generating the reference sequences of 262144 variants

1.fold_RNA.py: predicting RNA structure of each variant using RNAfold

2.process_RNA_structure.py: classifying variants into three structural groups (symmetric, asymmetric, others)

3.calculate_score.py: processing raw count to calculate global cleavage efficiency (gce), local cleavage efficiency (lce), accuracy (rc) scores

4.analysis_draw.py: analysing scores and generating figures

using_now.py: containing scripts to generate scores for different motifs and import tables

common_tool.py: containing common scripts to process RNA structures
