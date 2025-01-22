UCE analysis
1/13/23
Isabel novick
# Adapted from Brent Faircloth's UCE analysis using phyluce and illumiprocessor


### PATH TO LOGGING IN TO SCC###
ssh inovick@scc1.bu.edu
cd /projectnb/mullenl/novick

# print working directory:
[inovick@scc1 uce_analysis]$ pwd
/projectnb/mullenl/novick/uce_analysis

# clone github into working directory:
git clone https://github.com/faircloth-lab/phyluce/

# Make everything executable (chmod)
chmod +x *.py
chmod +x *.ini
chmod +x *.toml


# Here is the website for the general pipeline analysis: https://phyluce.readthedocs.io/en/latest/daily-use/daily-use-3-uce-processing.html#uce-processing
# Here is the website for the UCE tutorial: https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html
# Instructions for setting up conda environment are here: /projectnb/mullenl/software/phyluce-1.7.1/phyluce-1.7.1/README.txt
# Phyluce and illumiprocessor are found here: /projectnb/mullenl/software/


# Activate conda environment
module load miniconda/4.12.0 #you can do this everywhere/anywhere
conda activate /projectnb/mullenl/software/phyluce-1.7.1/phyluce-1.7.1

# But what worked the best so far is just
conda activate

# To exit conda environment:
conda deactivate

# For a 4 core interactive session:
qrsh -pe omp 4

# Do not do this just on its own, or you won't be in the conda environment anymore! Just add it to the end of the job if you submit a jobs

# To use however many cores is required without having to change the qsub file:
illumiprocessor --cores $NSLOTS

# Find the environment file for phyluce and what the exact name is for the environment. Then move into the folder/directory where the envrionment is, then source activate the environment name.

# To make .config files:
touch filename.conf


#### Trimming using Illumiprocessor
# Example of config file:
[adapters]
i7:CTGTCTCTTATACACATCTCCGAGCCCACGAGAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:CTGTCTCTTATACACATCTGACGCTGCCGACGA*GTGTAGATCTCGGTGGTCGCCGTATCATT

[tag sequences]
i7-N735:CAAGGTGA
i5-N535:CATGACAC

[tag map]
9269-JK-6_035:i7-N735,i5-N535

[names]
phereoeca_spp_22_08:9269-JK-6_0

# i7 is for R1 file, forward 5' to 3', i5 is for R2 file, is reverse complement of given i5 index
# To search for index in file itself, type "less 'filename'", then /CTGTCTCTTATACACATCTCCGAGCCCACGAGAC for R1 (i7) and /CTGTCTCTTATACACATCTGACGCTGCCGACGA for R2 (i5), and the 8 bases before that adapter sequence are the index

# Config file for 9269-JK-6_001, 9269-JK-6_002, 9269-JK-6_003, 9269-JK-6_004:
# N701 -> i7 for 9269-JK-6_001
# N501 -> i5 for 9269-JK-6_001
# N702 -> i7 for 9269-JK-6_002
# N502 -> i5 for 9269-JK-6_002
# N703 -> i7 for 9269-JK-6_003
# N503 -> i5 for 9269-JK-6_003
# N704 -> i7 for 9269-JK-6_004
# N504 -> i5 for 9269-JK-6_004

[adapters]
i7:CTGTCTCTTATACACATCTCCGAGCCCACGAGAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:CTGTCTCTTATACACATCTGACGCTGCCGACGA*GTGTAGATCTCGGTGGTCGCCGTATCATT

[tag sequences]
i7-N701:TGAGTCAG
i5-N501:GACGTGAC
i7-N702:GAATGCTC
i5-N502:GGAGCTGC
i7-N703:GAATATCC
i5-N503:GCTGCATG
i7-N704:CTTATGAA
i5-N504:GCAATCGT

[tag map]
9269-JK-6_001_S1_L005:i7-N701,i5-N501
9269-JK-6_002_S1_L005:i7-N702,i5-N502
9269-JK-6_003_S1_L005:i7-N703,i5-N503
9269-JK-6_004_S1_L005:i7-N704,i5-N504

[names]
9269-JK-6_001_S1_L005:tinea_murariella
9269-JK-6_002_S1_L005:tinea_apicimaculella_21_08
9269-JK-6_003_S1_L005:tinea_occidentalis
9269-JK-6_004_S1_L005:tinea_apicimaculella_21_12

# Command to run:
illumiprocessor --input raw_fastq --output cleaned_data2 --no-merge --config test.conf --r1-pattern "{}_(R1|READ1|Read1|read1)_\\d+.fastq" --r2-pattern "{}_(R2|READ2|Read2|read2)_\\d+.fastq"

# All trimming stuff found in
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/trimming_dir

# Corrupted data files are:
9269-JK-6_008_S1_L005_R1_001.fastq.gz
9269-JK-6_008_S1_L005_R2_001.fastq.gz
9269-JK-6_018_S1_L005_R1_001.fastq.gz
9269-JK-6_020_S1_L005_R1_001.fastq.gz
9269-JK-6_029_S1_L005_R2_001.fastq.gz
9269-JK-6_034_S1_L005_R1_001.fastq.gz

# Now downloading from Jen's new data sharing platform directly to the SCC, not using cyberduck. Made new directory found in:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/trimming_dir/2_corrupted-uce-data

# New directory for trimming corrupted files found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/trimming_dir/test_corrupted_files_dir

# Testing with just 9269-JK-6_008_S1_L005_R1_001.fastq.gz files since both R1 and R2 are there to see if the naming error was coming from not having the paired reads together

# Old corrupted files (that are better now) with paired end mate reads in
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/trimming_dir/test_corrupted_files_dir/2_corrupted-uce-data

# This worked. I ran them one at a time, with their complementary R1 or R2 file, and moved the output files to this directory:

/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/trimming_dir/trimmed_fastq

# Next I will make a .conf file for assembly (assmebly.conf)
# For assembly, we are running as an mpi job because it's so huge. Keep getting this error:

File "/projectnb/mullenl/software/phyluce-1.7.1/phyluce-1.7.1/lib/python3.6/os.py", line 220, in makedirs
   mkdir(name, mode)
FileExistsError: [Errno 17] File exists: '/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/assembly_dir/spades-assemblies'

# Ended up doing an array job

# Assembly output is found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array

# Contigs are found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array/contigs

# For assembly QC
# Run this script against all directories of reads

for i in spades-assemblies-array/contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

#Output:
# samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
10_monopis_marginistrigella.contigs.fasta,42287,15558884,367.93539385626786,1.4473031257803541,56,5212,268.0,1745
11_monopis_monachella_21_26.contigs.fasta,114705,31925255,278.3248768580271,0.7128549131752976,54,9874,250.0,1987
12_monopis_monachella_21_27.contigs.fasta,41602,14788780,355.482428729388,1.4671219759525505,56,4749,268.0,1776
13_monopis_pavlovski.contigs.fasta,176087,56709783,322.0554782579066,0.7369575286368376,52,32806,259.0,6074
14_morophaga_bucephala.contigs.fasta,45663,19761117,432.7599369292425,1.6159709094554247,54,9489,300.0,2965
15_niditinea_fuscella_22_09.contigs.fasta,27481,10978685,399.5009279138314,2.1392051658768314,56,8755,276.0,1625
16_niditinea_fuscella_22_10.contigs.fasta,23380,9721875,415.8201454234388,2.494134984052678,56,9215,278.0,1630
17_opogona_nipponica_21_31.contigs.fasta,38050,14511697,381.3849408672799,1.7891225994834992,56,10675,271.0,1943
18_opogona_nipponica_21_32.contigs.fasta,63067,22888339,362.9210046458528,1.1789003039475,56,6937,269.0,2419
19_phereoeca_spp_22_07.contigs.fasta,196198,59264656,302.0655460300309,0.5096384155444537,51,13607,260.0,2863
1_acrolophus_spp.contigs.fasta,37128,14758363,397.49954212454213,1.7612335994707067,52,10103,280.0,2088
20_phereoeca_spp_22_08.contigs.fasta,115585,36585813,316.5273435134317,0.7755206740311746,54,15612,263.0,2565
21_scardiella_approximatella.contigs.fasta,117508,35339106,300.7378731660823,0.8940205558163807,53,11052,255.0,2886
22_tinea_apicimaculella_21_08.contigs.fasta,153196,59792115,390.2981474712133,0.8782474454431511,56,10841,292.0,6079
23_tinea_apicimaculella_21_12.contigs.fasta,47036,18228690,387.5476230972021,1.4981214231784472,56,10296,275.0,2462
24_tinea_carnariella.contigs.fasta,173684,56868414,327.42459869648326,0.6094837021011758,56,9705,268.0,3520
25_tinea_murariella.contigs.fasta,68535,23499581,342.8843802436711,1.0342203349722896,56,9110,270.0,1943
26_tinea_occidentalis.contigs.fasta,165858,52620426,317.2619107911587,0.5823426463349778,52,9504,270.0,2598
27_tinea_pellionella_22_01.contigs.fasta,33798,14324712,423.8331262204864,1.8973928790527899,56,9149,291.0,2043
28_tinea_pellionella_22_02.contigs.fasta,35423,14297228,403.61426192022134,1.739929760513686,56,7840,283.0,1768
29_tinea_pellionella_22_03.contigs.fasta,28200,11783555,417.85656028368794,2.2270924478982606,56,9180,281.0,1898
2_doleromorpha_porphyria_22_20.contigs.fasta,38628,13915276,360.23806565185873,1.704968421118615,54,6100,266.0,1895
30_tineola_bisselliella_21_02.contigs.fasta,239624,97205283,405.6575426501519,0.6554670840309227,51,15843,324.0,10431
31_tineola_bisselliella_22_17.contigs.fasta,23354,10266298,439.5948445662413,2.678979509835316,56,6456,281.0,1918
32_tineola_bisselliella_22_24.contigs.fasta,24383,11082448,454.51535906164133,2.8020566069292636,56,8901,284.0,2180
33_trichophaga_tapetzella.contigs.fasta,226042,65391532,289.2893002185435,0.38151087451659726,53,15398,258.0,1841
34_xystrologa_wielgusi_21_37.contigs.fasta,32834,15480556,471.47944204178594,2.2161423372130584,56,8308,311.0,2833
35_xystrologa_wielgusi_21_38.contigs.fasta,65410,25831043,394.9096927075371,1.2446584348921668,56,12486,277.0,3186
3_doleromorpha_porphyria_22_24.contigs.fasta,36138,13146294,363.780342022248,1.789597708791611,50,6783,268.0,1856
4_erechthias_zebrina_22_16.contigs.fasta,51591,18766407,363.7535035180555,1.4424374906451252,53,7520,266.0,2306
5_erechthias_zebrina_22_26.contigs.fasta,127547,37511324,294.0980501305401,0.7401283106869104,54,15368,253.0,2725
6_hybroma_servulella.contigs.fasta,218196,66330500,303.99503198958735,0.6629031323217219,54,16655,256.0,4999
7_monopis_crocicapitella_21_22.contigs.fasta,34252,13854033,404.4736949667173,1.8860325318690685,56,13457,282.0,1927
8_monopis_crocicapitella_22_04.contigs.fasta,105839,32950561,311.3272139759446,0.8814095669989986,52,15509,257.0,2654
9_monopis_laevigella.contigs.fasta,167372,51117247,305.4109827211242,0.7071148668955611,55,14374,258.0,3465



# Match contigs to probes
# Get summary stats on exploded fastas
# samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
acrolophus-spp.unaligned.fasta,797,1056651,1325.7854454203261,15.816063759704209,150,3019,1349.0,630
doleromorpha-porphyria-22-20.unaligned.fasta,797,1076931,1351.2308657465496,16.147094810934608,123,3938,1374.0,629
doleromorpha-porphyria-22-24.unaligned.fasta,824,1104621,1340.5594660194174,16.118858703845994,230,3994,1360.0,650
erechthias-zebrina-22-16.unaligned.fasta,407,508971,1250.5429975429975,26.29509235028689,118,3255,1309.0,281
erechthias-zebrina-22-26.unaligned.fasta,679,841147,1238.80265095729,19.6804370896675,104,3554,1246.0,464
hybroma-servulella.unaligned.fasta,774,1076583,1390.9341085271317,19.1064518588621,56,5720,1386.0,620
monopis-crocicapitella-21-22.unaligned.fasta,718,996667,1388.1155988857938,18.52226358125989,111,4641,1411.5,576
monopis-crocicapitella-22-04.unaligned.fasta,708,955898,1350.138418079096,18.604637741151624,63,4392,1381.5,547
monopis-laevigella.unaligned.fasta,725,1007443,1389.576551724138,18.071914669762208,169,3251,1422.0,585
monopis-marginistrigella.unaligned.fasta,698,1060445,1519.2621776504297,18.376268427281275,208,3104,1574.5,604
monopis-monachella-21-26.unaligned.fasta,652,805135,1234.8696319018404,17.72084671236989,65,3133,1244.0,460
monopis-monachella-21-27.unaligned.fasta,637,778051,1221.4301412872842,17.375620051787973,125,3231,1248.0,442
monopis-pavlovski.unaligned.fasta,635,731535,1152.023622047244,17.217697541650473,73,2519,1187.0,416
morophaga-bucephala.unaligned.fasta,772,1205368,1561.3575129533679,14.978068600634186,99,3270,1619.0,703
niditinea-fuscella-22-09.unaligned.fasta,678,915275,1349.9631268436578,18.50905487450637,139,3419,1377.0,528
niditinea-fuscella-22-10.unaligned.fasta,700,961041,1372.9157142857143,18.443340138855117,163,3180,1386.0,553
opogona-nipponica-21-31.unaligned.fasta,792,1209745,1527.455808080808,17.668237227665774,125,3957,1563.0,688
opogona-nipponica-21-32.unaligned.fasta,776,1144314,1474.631443298969,18.699540486894062,186,3765,1514.0,646
phereoeca-spp-22-07.unaligned.fasta,703,972264,1383.0213371266002,18.94710831018974,111,4667,1364.0,559
phereoeca-spp-22-08.unaligned.fasta,756,1104640,1461.164021164021,20.452120040853593,102,6963,1475.0,634
scardiella-approximatella.unaligned.fasta,727,992996,1365.881705639615,19.85262591040319,86,5935,1367.0,564
tinea-apicimaculella-21-08.unaligned.fasta,663,1018933,1536.8521870286577,25.915879029435565,111,4861,1545.0,525
tinea-apicimaculella-21-12.unaligned.fasta,678,910516,1342.94395280236,21.15582134061362,68,3709,1354.5,498
tinea-carnariella.unaligned.fasta,786,1093298,1390.9643765903309,19.67824394390615,108,5802,1392.5,614
tinea-murariella.unaligned.fasta,672,1040953,1549.0372023809523,18.332974465398763,130,3770,1602.0,586
tinea-occidentalis.unaligned.fasta,513,678184,1321.9961013645225,31.1151827508525,202,4453,1287.0,343
tinea-pellionella-22-01.unaligned.fasta,474,622514,1313.3206751054852,22.531590857823332,148,3226,1326.5,361
tinea-pellionella-22-02.unaligned.fasta,708,970680,1371.0169491525423,20.096919022321842,111,6955,1376.0,561
tinea-pellionella-22-03.unaligned.fasta,287,351560,1224.9477351916375,32.31595352938979,242,3656,1258.0,188
tineola-bisselliella-21-02.unaligned.fasta,592,1023238,1728.4425675675675,34.85517125258941,77,7771,1694.0,479
tineola-bisselliella-22-17.unaligned.fasta,660,998832,1513.3818181818183,22.99025735172184,194,4545,1559.0,540
tineola-bisselliella-22-24.unaligned.fasta,666,1045298,1569.5165165165165,24.604193544807202,111,4517,1590.0,542
trichophaga-tapetzella.unaligned.fasta,653,711335,1089.3338437978562,14.80922736447562,59,2392,1091.0,396
xystrologa-wielgusi-21-37.unaligned.fasta,830,1377699,1659.878313253012,16.619960291123935,185,3614,1708.5,756
xystrologa-wielgusi-21-38.unaligned.fasta,809,1290906,1595.68108776267,20.792273258679195,228,12486,1637.0,727

# Aligning UCE loci
# Edge trimming found in
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/taxon_sets/all

# Making an edge trimming qsub in this directory:
edge_trimming.qsub

# Almost done, got an error when creating the final data matrices for 95p and 75p. The error was:

Traceback (most recent call last):
File "/projectnb/mullenl/software/phyluce-1.7.1/phyluce-1.7.1/lib/python3.6/logging/__init__.py", line 996, in emit
stream.write(msg)
UnicodeEncodeError: 'ascii' codec cant encode character '\u2265' in position 132: ordinal not in range(128)

# The fix is to add
export LC_ALL=en_US.UTF-8
# To the qsub file before calling python and submitting the job. Let's see if it works

# It works!

# Now to use the phylip file in raxml. The phylip file for 75p is found here: /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/taxon_sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml

# The 95p is found in the 95p-raxml directory
# To use raxml: use qsub found in raxml directory
rojectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/raxml/raxml.qsub

# To replace uce- with uce_ in the nexus file, use this command:
sed -i 's/uce-/uce_/g' name-of-nexus-file.nexus


###### Doing one more time to exclude tinea_apicimaculella_21_08 (poor quality)
# Contigs found in
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/match_contigs_probes_dir/2_uce-search-results/probe.matches.sqlite

# For the no_apic_95 qsub:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/no_apic_mafft-nexus-internal-trimmed-gblocks-clean-95p-raxml/no_apic_95p_raxml.phylip

######### 7/5/23: Now have to harvest uce loci from genomes (T. pellionella, M. laevigella, T. trinotella, T. semifulvella)
# Using this tutorial: https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-3.html#tutorial-iii-harvesting-uce-loci-from-genomes
# The genomes can be found here: /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes

##### 11/8/23 Harvesting uce loci from more genomes!
# adela_reaumurella_genome: adela_reaumurella_genome.2bit found here: /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/adela_reaumurella_genome/data/GCA_009867175.2
# cydia_amplana_genome: cydia_amplana_genome.2bit found here: /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/cydia_amplana_genome/data/GCA_948474715.1

# Genomes are in .fna file format. I just changed the suffix to .fa (fasta) format and hope that works. Now I need to change it from a fasta to a 2bit
# Used this command to download program:
rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./
from http://hgdownload.soe.ucsc.edu/downloads.html 
# Source_downloads that will allow conversion from .fa to 2bit. Not sure if it will work
# This downloaded a ton of stuff into this directory: /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/fa_to_2bit_dir/fa_to_2bit

# I'm going to copy the adela genome and the cydia genome into this 2bit directory and use the following command
# Converting to 2bit:
faToTwoBit genome.fa genome.2bit

# The 2bit genome is now called tin_pell_genome.2bit in tin_pell_genome directory
# The 2bit genomes for adela and cydia will also be in their respective directories, as well as 2 new directories that only contain the 2bit genomes:
cydia 2bit genome: /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/cydia_amp_genome
adela 2bit genome: /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/adela_rea_genome

# This is so I can use the align_probes_to_genome.qsub that's already made in the genomes directory that will hopefully work with the scaffoldlist

# Now made a qsub called align_probes_to_genome.qsub, ran, created aligned-probes-genome-lastz where cleaned alignments are located
# Got an error when running with outgroups "this isn't a file"
# This was the fix:
# 1) first, you need to give the exact organism name you used for the 2bit file, in --scaffoldlist option, i.e. :
# it shall be:
--scaffoldlist cydia_amplana_genome adela_reaumurella_genome
rather than:
--scaffoldlist cydia_amp_genome adela_rea_genome

# since your 2bit file names respectively are: 
--scaffoldlist cydia_amp_genome.2bit and adela_rea_genome.2bit

# The program will look for the 2bit file name that matches with what is told in --scaffoldlist.

# 2) further more, simply made above change will not work, since your genome file is put in the subfolder, not in the base folder you give in
--genome-base-path /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/ \

# so you need to either move them or create a symbolic link for each genome file under that base folder. I did the latter. see below: (you use 'ln -s' command to link a file, please refer to the command help for more detail)
lrwxrwxrwx 1 inovick mullenl 42 Nov 13 15:20 cydia_amplana_genome.2bit -> cydia_amp_genome/cydia_amplana_genome.2bit
lrwxrwxrwx 1 inovick mullenl 46 Nov 13 15:05 adela_reaumurella_genome.2bit -> adela_rea_genome/adela_reaumurella_genome.2bit

# 3) even with this, it seems to me you would still need to add one more option '--no-dir' to indicate that your actual genome seq. data is not in the same folder as the 2bit files is located.

# So after made the above changes, I finally run the following command:
# This is the new qsub
phyluce_probe_run_multiple_lastzs_sqlite \
--db align_probes_to_genome_outgroups.sqlite \
--output aligned-probes-genome-outgroups-lastz \
--probefile Lepidoptera-UCE-1.3K-v1.fasta \
--scaffoldlist cydia_amplana_genome adela_reaumurella_genome \
--genome-base-path /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/ \
--no-dir \
--cores 16

# Now made a config file called scaffolds.conf that contains the path to the genomes/scaffolds that are Used
# Made qsub to extract FASTA data from each genome for each UCE locus called extract_sequences.qsub. Output found in extract_sequences.qlog. Really small amount of UCEs pulled compared to the tutorial, probably because they weren't specifically targeted in sequencing like other samples

# New config file for outgroups is called scaffolds_2.conf and qsub to extract FASTA data is called outgroups_extract_sequences.qsub
# It is not working, getting an error that file does not exist (the lastz cleaned file)

# Naming scheme is what is giving issues. Name of directory needs to match name of genome, so I am going through, renaming everything, reorganizing, and starting from beginning
# After renaming, aligning probes to genome worked. Extracting is still giving issues
# It worked after I adjusted the path name to {}_genome

# Use module load miniconda 23.1.0 (newer version)
# Next step is to symlink the following (found in /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/extract-sequences-fasta)
-rw-r--r-- 1 inovick mullenl 74K Jul 17 10:22 mon_laev_genome.fasta
-rw-r--r-- 1 inovick mullenl 67K Jul 17 10:22 tin_pell_genome.fasta
-rw-r--r-- 1 inovick mullenl 71K Jul 17 10:22 tin_semi_genome.fasta
-rw-r--r-- 1 inovick mullenl 73K Jul 17 10:22 tin_trino_genome.fasta

# into /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array/contigs

# This is the command i used to try and create the symlink:
ln -s /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/extract-sequences-fasta/mon_laev_genome.fasta ~/mon_laev_genome.fasta

# But permission is denied

# Do This:
ln -s /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/extract-sequences-fasta/mon_laev_genome.fasta /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array/contigs/mon_laev_genome.fasta

# mon_laev_genome.fasta got deleted, reran extract_sequences.qsub, now it can be found in
/projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/extract-sequences-fasta-2

# To symlink, I will use
ln -s /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/extract-sequences-fasta-2/mon_laev_genome.fasta /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array/contigs/mon_laev_genome.fasta

# For cydia_amp: 
ln -s /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/outgroups_extract-sequences-fasta/cydia_amp.fasta /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array/contigs/cydia_amp.fasta

# For adela_rea: 
ln -s /projectnb/mullenl/novick/uce_analysis/uce_genomes/genomes/outgroups_extract-sequences-fasta/adela_rea.fasta /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array/contigs/adela_rea.fasta

# After qc:
# samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb:
acrolophus_spp.contigs.fasta,37128,14758363,397.49954212454213,1.7612335994707067,52,10103,280.0,2088
doleromorpha_porphyria_22_20.contigs.fasta,38628,13915276,360.23806565185873,1.704968421118615,54,6100,266.0,1895
doleromorpha_porphyria_22_24.contigs.fasta,36138,13146294,363.780342022248,1.789597708791611,50,6783,268.0,1856
erechthias_zebrina_22_16.contigs.fasta,51591,18766407,363.7535035180555,1.4424374906451252,53,7520,266.0,2306
erechthias_zebrina_22_26.contigs.fasta,127547,37511324,294.0980501305401,0.7401283106869104,54,15368,253.0,2725
hybroma_servulella.contigs.fasta,218196,66330500,303.99503198958735,0.6629031323217219,54,16655,256.0,4999
mon_laev_genome.fasta,61,65894,1080.2295081967213,12.44854751425033,590,1190,1109.0,53
monopis_crocicapitella_21_22.contigs.fasta,34252,13854033,404.4736949667173,1.8860325318690685,56,13457,282.0,1927
monopis_crocicapitella_22_04.contigs.fasta,105839,32950561,311.3272139759446,0.8814095669989986,52,15509,257.0,2654
monopis_laevigella.contigs.fasta,167372,51117247,305.4109827211242,0.7071148668955611,55,14374,258.0,3465
monopis_marginistrigella.contigs.fasta,42287,15558884,367.93539385626786,1.4473031257803541,56,5212,268.0,1745
monopis_monachella_21_26.contigs.fasta,114705,31925255,278.3248768580271,0.7128549131752976,54,9874,250.0,1987
monopis_monachella_21_27.contigs.fasta,41602,14788780,355.482428729388,1.4671219759525505,56,4749,268.0,1776
monopis_pavlovski.contigs.fasta,176087,56709783,322.0554782579066,0.7369575286368376,52,32806,259.0,6074
morophaga_bucephala.contigs.fasta,45663,19761117,432.7599369292425,1.6159709094554247,54,9489,300.0,2965
niditinea_fuscella_22_09.contigs.fasta,27481,10978685,399.5009279138314,2.1392051658768314,56,8755,276.0,1625
niditinea_fuscella_22_10.contigs.fasta,23380,9721875,415.8201454234388,2.494134984052678,56,9215,278.0,1630
opogona_nipponica_21_31.contigs.fasta,38050,14511697,381.3849408672799,1.7891225994834992,56,10675,271.0,1943
opogona_nipponica_21_32.contigs.fasta,63067,22888339,362.9210046458528,1.1789003039475,56,6937,269.0,2419
phereoeca_spp_22_07.contigs.fasta,196198,59264656,302.0655460300309,0.5096384155444537,51,13607,260.0,2863
phereoeca_spp_22_08.contigs.fasta,115585,36585813,316.5273435134317,0.7755206740311746,54,15612,263.0,2565
scardiella_approximatella.contigs.fasta,117508,35339106,300.7378731660823,0.8940205558163807,53,11052,255.0,2886
tin_pell_genome.fasta,54,59911,1109.462962962963,7.828691383567851,852,1185,1121.0,52
tin_semi_genome.fasta,59,63190,1071.0169491525423,13.249118968214914,759,1193,1108.0,49
tin_trino_genome.fasta,60,65637,1093.95,7.7549118751982205,854,1199,1108.0,57
tinea_apicimaculella_21_12.contigs.fasta,47036,18228690,387.5476230972021,1.4981214231784472,56,10296,275.0,2462
tinea_carnariella.contigs.fasta,173684,56868414,327.42459869648326,0.6094837021011758,56,9705,268.0,3520
tinea_murariella.contigs.fasta,68535,23499581,342.8843802436711,1.0342203349722896,56,9110,270.0,1943
tinea_occidentalis.contigs.fasta,165858,52620426,317.2619107911587,0.5823426463349778,52,9504,270.0,2598
tinea_pellionella_22_01.contigs.fasta,33798,14324712,423.8331262204864,1.8973928790527899,56,9149,291.0,2043
tinea_pellionella_22_02.contigs.fasta,35423,14297228,403.61426192022134,1.739929760513686,56,7840,283.0,1768
tinea_pellionella_22_03.contigs.fasta,28200,11783555,417.85656028368794,2.2270924478982606,56,9180,281.0,1898
tineola_bisselliella_21_02.contigs.fasta,239624,97205283,405.6575426501519,0.6554670840309227,51,15843,324.0,10431
tineola_bisselliella_22_17.contigs.fasta,23354,10266298,439.5948445662413,2.678979509835316,56,6456,281.0,1918
tineola_bisselliella_22_24.contigs.fasta,24383,11082448,454.51535906164133,2.8020566069292636,56,8901,284.0,2180
trichophaga_tapetzella.contigs.fasta,226042,65391532,289.2893002185435,0.38151087451659726,53,15398,258.0,1841
xystrologa_wielgusi_21_37.contigs.fasta,32834,15480556,471.47944204178594,2.2161423372130584,56,8308,311.0,2833
xystrologa_wielgusi_21_38.contigs.fasta,65410,25831043,394.9096927075371,1.2446584348921668,56,12486,277.0,3186


# Outgroups: I think what I did before is start at "Assembly QC"
# We can get a sense of how well the assembly worked by running the following from the top of our working directory:

for i in spades-assemblies/contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

# For the run including outgroups (this is the final run containing everything, all NCBI genomes, and getting rid of bad sample:
# samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
adela_rea.fasta,48,50255,1046.9791666666667,23.478005581490848,340,1170,1110.5,39
cydia_amp.fasta,174,194651,1118.683908045977,20.769386236853236,285,3661,1120.0,155
doleromorpha_porphyria_22_20.contigs.fasta,38628,13915276,360.23806565185873,1.704968421118615,54,6100,266.0,1895
doleromorpha_porphyria_22_24.contigs.fasta,36138,13146294,363.780342022248,1.789597708791611,50,6783,268.0,1856
erechthias_zebrina_22_16.contigs.fasta,51591,18766407,363.7535035180555,1.4424374906451252,53,7520,266.0,2306
erechthias_zebrina_22_26.contigs.fasta,127547,37511324,294.0980501305401,0.7401283106869104,54,15368,253.0,2725
hybroma_servulella.contigs.fasta,218196,66330500,303.99503198958735,0.6629031323217219,54,16655,256.0,4999
mon_laev_genome.fasta,61,65894,1080.2295081967213,12.44854751425033,590,1190,1109.0,53
monopis_crocicapitella_21_22.contigs.fasta,34252,13854033,404.4736949667173,1.8860325318690685,56,13457,282.0,1927
monopis_crocicapitella_22_04.contigs.fasta,105839,32950561,311.3272139759446,0.8814095669989986,52,15509,257.0,2654
monopis_laevigella.contigs.fasta,167372,51117247,305.4109827211242,0.7071148668955611,55,14374,258.0,3465
monopis_marginistrigella.contigs.fasta,42287,15558884,367.93539385626786,1.4473031257803541,56,5212,268.0,1745
monopis_monachella_21_26.contigs.fasta,114705,31925255,278.3248768580271,0.7128549131752976,54,9874,250.0,1987
monopis_monachella_21_27.contigs.fasta,41602,14788780,355.482428729388,1.4671219759525505,56,4749,268.0,1776
monopis_pavlovski.contigs.fasta,176087,56709783,322.0554782579066,0.7369575286368376,52,32806,259.0,6074
morophaga_bucephala.contigs.fasta,45663,19761117,432.7599369292425,1.6159709094554247,54,9489,300.0,2965
niditinea_fuscella_22_09.contigs.fasta,27481,10978685,399.5009279138314,2.1392051658768314,56,8755,276.0,1625
niditinea_fuscella_22_10.contigs.fasta,23380,9721875,415.8201454234388,2.494134984052678,56,9215,278.0,1630
opogona_nipponica_21_31.contigs.fasta,38050,14511697,381.3849408672799,1.7891225994834992,56,10675,271.0,1943
opogona_nipponica_21_32.contigs.fasta,63067,22888339,362.9210046458528,1.1789003039475,56,6937,269.0,2419
phereoeca_spp_22_07.contigs.fasta,196198,59264656,302.0655460300309,0.5096384155444537,51,13607,260.0,2863
phereoeca_spp_22_08.contigs.fasta,115585,36585813,316.5273435134317,0.7755206740311746,54,15612,263.0,2565
scardiella_approximatella.contigs.fasta,117508,35339106,300.7378731660823,0.8940205558163807,53,11052,255.0,2886
tinea_apicimaculella_21_12.contigs.fasta,47036,18228690,387.5476230972021,1.4981214231784472,56,10296,275.0,2462
tinea_carnariella.contigs.fasta,173684,56868414,327.42459869648326,0.6094837021011758,56,9705,268.0,3520
tinea_murariella.contigs.fasta,68535,23499581,342.8843802436711,1.0342203349722896,56,9110,270.0,1943
tinea_occidentalis.contigs.fasta,165858,52620426,317.2619107911587,0.5823426463349778,52,9504,270.0,2598
tinea_pellionella_22_01.contigs.fasta,33798,14324712,423.8331262204864,1.8973928790527899,56,9149,291.0,2043
tinea_pellionella_22_02.contigs.fasta,35423,14297228,403.61426192022134,1.739929760513686,56,7840,283.0,1768
tinea_pellionella_22_03.contigs.fasta,28200,11783555,417.85656028368794,2.2270924478982606,56,9180,281.0,1898
tineola_bisselliella_21_02.contigs.fasta,239624,97205283,405.6575426501519,0.6554670840309227,51,15843,324.0,10431
tineola_bisselliella_22_17.contigs.fasta,23354,10266298,439.5948445662413,2.678979509835316,56,6456,281.0,1918
tineola_bisselliella_22_24.contigs.fasta,24383,11082448,454.51535906164133,2.8020566069292636,56,8901,284.0,2180
tin_pell_genome.fasta,54,59911,1109.462962962963,7.828691383567851,852,1185,1121.0,52
tin_semi_genome.fasta,59,63190,1071.0169491525423,13.249118968214914,759,1193,1108.0,49
tin_trino_genome.fasta,60,65637,1093.95,7.7549118751982205,854,1199,1108.0,57
trichophaga_tapetzella.contigs.fasta,226042,65391532,289.2893002185435,0.38151087451659726,53,15398,258.0,1841
xystrologa_wielgusi_21_37.contigs.fasta,32834,15480556,471.47944204178594,2.2161423372130584,56,8308,311.0,2833
xystrologa_wielgusi_21_38.contigs.fasta,65410,25831043,394.9096927075371,1.2446584348921668,56,12486,277.0,3186


# Now running phyluce_assembly_match_contigs_to_probes program
# It's the match_contigs.qsub located in /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/match_contigs_probes_dir
# New output is in outgroups_uce-search-results

adela_rea: 45 (93.75%) uniques of 48 contigs, 0 dupe probe matches, 3 UCE loci removed for matching multiple
cydia_amp: 150 (86.21%) uniques of 174 contigs, 0 dupe probe matches, 16 UCE loci removed for matching multiple

# Total info found in phyluce_assembly_match_contigs_to_probes.log, current run directory is /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/match_contigs_probes_dir/outgroups_uce-search-results

# Everything at this point is in
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all

# But all final run stuff has prefix 3

# All outgroups related stuff will have outgroup in the name

# Now doing extracting uce loci and making taxon set conf file called outgroups_taxon_set.conf
# Used in qsub called outgroups_data_matrix.qsub to make the data matrix configuration file

# Making new qsub for extracting FASTA data that correspond to the loci in all-taxa-incomplete.conf called get_fasta_data.qsub

# After exploding fastas, here are summary stats on all fastas including outgroups:
samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
acrolophus-spp.unaligned.fasta,797,1056651,1325.7854454203261,15.816063759704209,150,3019,1349.0,630
adela-rea.unaligned.fasta,45,46628,1036.1777777777777,25.471439519605045,340,1170,1105.0,36
cydia-amp.unaligned.fasta,150,167642,1117.6133333333332,24.079186639750034,285,3661,1120.0,131
doleromorpha-porphyria-22-20.unaligned.fasta,797,1076931,1351.2308657465496,16.147094810934608,123,3938,1374.0,629
doleromorpha-porphyria-22-24.unaligned.fasta,824,1104621,1340.5594660194174,16.118858703845994,230,3994,1360.0,650
erechthias-zebrina-22-16.unaligned.fasta,407,508971,1250.5429975429975,26.29509235028689,118,3255,1309.0,281
erechthias-zebrina-22-26.unaligned.fasta,679,841147,1238.80265095729,19.6804370896675,104,3554,1246.0,464
hybroma-servulella.unaligned.fasta,774,1076583,1390.9341085271317,19.1064518588621,56,5720,1386.0,620
mon-laev-genome.unaligned.fasta,52,55754,1072.1923076923076,14.282805071485754,590,1190,1102.5,44
monopis-crocicapitella-21-22.unaligned.fasta,718,996667,1388.1155988857938,18.52226358125989,111,4641,1411.5,576
monopis-crocicapitella-22-04.unaligned.fasta,708,955898,1350.138418079096,18.604637741151624,63,4392,1381.5,547
monopis-laevigella.unaligned.fasta,725,1007443,1389.576551724138,18.071914669762208,169,3251,1422.0,585
monopis-marginistrigella.unaligned.fasta,698,1060445,1519.2621776504297,18.376268427281275,208,3104,1574.5,604
monopis-monachella-21-26.unaligned.fasta,652,805135,1234.8696319018404,17.72084671236989,65,3133,1244.0,460
monopis-monachella-21-27.unaligned.fasta,637,778051,1221.4301412872842,17.375620051787973,125,3231,1248.0,442
monopis-pavlovski.unaligned.fasta,635,731535,1152.023622047244,17.217697541650473,73,2519,1187.0,416
morophaga-bucephala.unaligned.fasta,772,1205368,1561.3575129533679,14.978068600634186,99,3270,1619.0,703
niditinea-fuscella-22-09.unaligned.fasta,678,915275,1349.9631268436578,18.50905487450637,139,3419,1377.0,528
niditinea-fuscella-22-10.unaligned.fasta,700,961041,1372.9157142857143,18.443340138855117,163,3180,1386.0,553
opogona-nipponica-21-31.unaligned.fasta,792,1209745,1527.455808080808,17.668237227665774,125,3957,1563.0,688
opogona-nipponica-21-32.unaligned.fasta,776,1144314,1474.631443298969,18.699540486894062,186,3765,1514.0,646
phereoeca-spp-22-07.unaligned.fasta,703,972264,1383.0213371266002,18.94710831018974,111,4667,1364.0,559
phereoeca-spp-22-08.unaligned.fasta,756,1104640,1461.164021164021,20.452120040853593,102,6963,1475.0,634
scardiella-approximatella.unaligned.fasta,727,992996,1365.881705639615,19.85262591040319,86,5935,1367.0,564
tinea-apicimaculella-21-12.unaligned.fasta,678,910516,1342.94395280236,21.15582134061362,68,3709,1354.5,498
tinea-carnariella.unaligned.fasta,786,1093298,1390.9643765903309,19.67824394390615,108,5802,1392.5,614
tinea-murariella.unaligned.fasta,672,1040953,1549.0372023809523,18.332974465398763,130,3770,1602.0,586
tinea-occidentalis.unaligned.fasta,513,678184,1321.9961013645225,31.1151827508525,202,4453,1287.0,343
tinea-pellionella-22-01.unaligned.fasta,474,622514,1313.3206751054852,22.531590857823332,148,3226,1326.5,361
tinea-pellionella-22-02.unaligned.fasta,708,970680,1371.0169491525423,20.096919022321842,111,6955,1376.0,561
tinea-pellionella-22-03.unaligned.fasta,287,351560,1224.9477351916375,32.31595352938979,242,3656,1258.0,188
tineola-bisselliella-21-02.unaligned.fasta,592,1023238,1728.4425675675675,34.85517125258941,77,7771,1694.0,479
tineola-bisselliella-22-17.unaligned.fasta,660,998832,1513.3818181818183,22.99025735172184,194,4545,1559.0,540
tineola-bisselliella-22-24.unaligned.fasta,666,1045298,1569.5165165165165,24.604193544807202,111,4517,1590.0,542
tin-pell-genome.unaligned.fasta,43,47645,1108.0232558139535,9.298572055564675,852,1185,1121.0,41
tin-semi-genome.unaligned.fasta,52,55287,1063.2115384615386,14.665946863039164,759,1193,1098.5,42
tin-trino-genome.unaligned.fasta,50,54461,1089.22,8.933326276776581,854,1199,1102.5,47
trichophaga-tapetzella.unaligned.fasta,653,711335,1089.3338437978562,14.80922736447562,59,2392,1091.0,396
xystrologa-wielgusi-21-37.unaligned.fasta,830,1377699,1659.878313253012,16.619960291123935,185,3614,1708.5,756
xystrologa-wielgusi-21-38.unaligned.fasta,809,1290906,1595.68108776267,20.792273258679195,228,12486,1637.0,727

# Summary stats on alignments after edge trimming after running this including outgroups:
phyluce_align_get_align_summary_data \
    --alignments outgroups_mafft-nexus-edge-trimmed \
    --cores 12 \
    --log-path log
    
# Gives informative sites summary:
2023-11-29 09:30:46,598 - phyluce_align_get_align_summary_data - INFO - ------------------- Informative Sites summary -------------------
2023-11-29 09:30:46,598 - phyluce_align_get_align_summary_data - INFO - [Sites] loci:   1,059
2023-11-29 09:30:46,598 - phyluce_align_get_align_summary_data - INFO - [Sites] total:  248,404
2023-11-29 09:30:46,598 - phyluce_align_get_align_summary_data - INFO - [Sites] mean:   234.56
2023-11-29 09:30:46,598 - phyluce_align_get_align_summary_data - INFO - [Sites] 95% CI: 13.63
2023-11-29 09:30:46,598 - phyluce_align_get_align_summary_data - INFO - [Sites] min:    0
2023-11-29 09:30:46,598 - phyluce_align_get_align_summary_data - INFO - [Sites] max:    1,187

# After running gblocks trumming, running this for summary stats:
phyluce_align_get_align_summary_data \
    --alignments outgroups_mafft-nexus-internal-trimmed-gblocks \
    --cores 12 \
    --log-path log
    
# Data matrix completeness summary:
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]		745 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]		652 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]		548 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]		429 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]		312 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]		189 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]		75 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]		45 alignments
2023-11-29 10:21:15,741 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]		13 alignments
2023-11-29 10:21:15,742 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]		4 alignments

# Creating the 75p data matrix:
# Copied 128 alignments of 1172 total containing ≥ 0.75 proportion of taxa (n = 30)

# Creating concatenated data matrix in nexus format:
outgroups_mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml-nexus

# Summary data for alignments is in /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/log/phyluce_align_get_align_summary_data.log
# Phylip file found in
3_mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml/3_mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml.phylip


# After exploding fastas and matching to probes, here are the summary stats of all samples:

samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
acrolophus-spp.unaligned.fasta,797,1056651,1325.7854454203261,15.816063759704209,150,3019,1349.0,630
doleromorpha-porphyria-22-20.unaligned.fasta,797,1076931,1351.2308657465496,16.147094810934608,123,3938,1374.0,629
doleromorpha-porphyria-22-24.unaligned.fasta,824,1104621,1340.5594660194174,16.118858703845994,230,3994,1360.0,650
erechthias-zebrina-22-16.unaligned.fasta,407,508971,1250.5429975429975,26.29509235028689,118,3255,1309.0,281
erechthias-zebrina-22-26.unaligned.fasta,679,841147,1238.80265095729,19.6804370896675,104,3554,1246.0,464
hybroma-servulella.unaligned.fasta,774,1076583,1390.9341085271317,19.1064518588621,56,5720,1386.0,620
mon-laev-genome.unaligned.fasta,52,55754,1072.1923076923076,14.282805071485754,590,1190,1102.5,44
monopis-crocicapitella-21-22.unaligned.fasta,718,996667,1388.1155988857938,18.52226358125989,111,4641,1411.5,576
monopis-crocicapitella-22-04.unaligned.fasta,708,955898,1350.138418079096,18.604637741151624,63,4392,1381.5,547
monopis-laevigella.unaligned.fasta,725,1007443,1389.576551724138,18.071914669762208,169,3251,1422.0,585
monopis-marginistrigella.unaligned.fasta,698,1060445,1519.2621776504297,18.376268427281275,208,3104,1574.5,604
monopis-monachella-21-26.unaligned.fasta,652,805135,1234.8696319018404,17.72084671236989,65,3133,1244.0,460
monopis-monachella-21-27.unaligned.fasta,637,778051,1221.4301412872842,17.375620051787973,125,3231,1248.0,442
monopis-pavlovski.unaligned.fasta,635,731535,1152.023622047244,17.217697541650473,73,2519,1187.0,416
morophaga-bucephala.unaligned.fasta,772,1205368,1561.3575129533679,14.978068600634186,99,3270,1619.0,703
niditinea-fuscella-22-09.unaligned.fasta,678,915275,1349.9631268436578,18.50905487450637,139,3419,1377.0,528
niditinea-fuscella-22-10.unaligned.fasta,700,961041,1372.9157142857143,18.443340138855117,163,3180,1386.0,553
opogona-nipponica-21-31.unaligned.fasta,792,1209745,1527.455808080808,17.668237227665774,125,3957,1563.0,688
opogona-nipponica-21-32.unaligned.fasta,776,1144314,1474.631443298969,18.699540486894062,186,3765,1514.0,646
phereoeca-spp-22-07.unaligned.fasta,703,972264,1383.0213371266002,18.94710831018974,111,4667,1364.0,559
phereoeca-spp-22-08.unaligned.fasta,756,1104640,1461.164021164021,20.452120040853593,102,6963,1475.0,634
scardiella-approximatella.unaligned.fasta,727,992996,1365.881705639615,19.85262591040319,86,5935,1367.0,564
tinea-apicimaculella-21-12.unaligned.fasta,678,910516,1342.94395280236,21.15582134061362,68,3709,1354.5,498
tinea-carnariella.unaligned.fasta,786,1093298,1390.9643765903309,19.67824394390615,108,5802,1392.5,614
tinea-murariella.unaligned.fasta,672,1040953,1549.0372023809523,18.332974465398763,130,3770,1602.0,586
tinea-occidentalis.unaligned.fasta,513,678184,1321.9961013645225,31.1151827508525,202,4453,1287.0,343
tinea-pellionella-22-01.unaligned.fasta,474,622514,1313.3206751054852,22.531590857823332,148,3226,1326.5,361
tinea-pellionella-22-02.unaligned.fasta,708,970680,1371.0169491525423,20.096919022321842,111,6955,1376.0,561
tinea-pellionella-22-03.unaligned.fasta,287,351560,1224.9477351916375,32.31595352938979,242,3656,1258.0,188
tineola-bisselliella-21-02.unaligned.fasta,592,1023238,1728.4425675675675,34.85517125258941,77,7771,1694.0,479
tineola-bisselliella-22-17.unaligned.fasta,660,998832,1513.3818181818183,22.99025735172184,194,4545,1559.0,540
tineola-bisselliella-22-24.unaligned.fasta,666,1045298,1569.5165165165165,24.604193544807202,111,4517,1590.0,542
tin-pell-genome.unaligned.fasta,43,47645,1108.0232558139535,9.298572055564675,852,1185,1121.0,41
tin-semi-genome.unaligned.fasta,52,55287,1063.2115384615386,14.665946863039164,759,1193,1098.5,42
tin-trino-genome.unaligned.fasta,50,54461,1089.22,8.933326276776581,854,1199,1102.5,47
trichophaga-tapetzella.unaligned.fasta,653,711335,1089.3338437978562,14.80922736447562,59,2392,1091.0,396
xystrologa-wielgusi-21-37.unaligned.fasta,830,1377699,1659.878313253012,16.619960291123935,185,3614,1708.5,756
xystrologa-wielgusi-21-38.unaligned.fasta,809,1290906,1595.68108776267,20.792273258679195,228,12486,1637.0,727

# Raxml file for future use that created figtree max likelihood tree is RAxML_bipartitions.outgroups_raxml_75p
# This was from outgroups_mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml-phylip which also contains a charsets file
# 1172 loci, 128 alignments

# For subsetting loci:

# R script found in /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/subsetting_loci/phyloch_script.R
# Copied r script to correct folder, /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/outgroups_mafft-nexus-internal-trimmed-gblocks-clean-75p
# This is where I'm going to subset the loci
# Going to make a directory full of the top 75 parsimony informative sites, copy
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/outgroups_mafft-nexus-internal-trimmed-gblocks-clean-75p/top75names.txt
# into
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/subsetting_loci/n75_pis_loci

#### For BEAST
# Get SCC to install BEAST and use interactive Rstudio session on SCC
# Filter for top 75 most informative sites, create new file
# Nexus files for sites are in
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/3_mafft-nexus-internal-trimmed-gblocks-clean-75p

# Nexus file (analogous to phylip file that we use for raxml) is located here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/3_mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml-nexus

# Copying over the sites:
for i in $(cat /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/outgroups_mafft-nexus-internal-trimmed-gblocks-clean-75p/top75names.txt); \
	do cp /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/outgroups_mafft-nexus-internal-trimmed-gblocks-clean-75p/$i n75_pis_loci/; \
done

# 75 most informative sites found Here
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/subsetting_loci/n75_pis_loci

## Iqtree:

# Now going to concatenate a data matrix with charsets to use in iqtree using concat_data_matrix_75p.qsub
# Data matrix of 75 most parsimony informative sites found in 
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/subsetting_loci/concat_matrix_75pis_loci_phylip

# Phylip file is in this directory along with charsets file
# Modelfinder also takes nexus files, which I will use because it has the "sets" section. Each of my loci is its own partition


# For modelfinder.qsub:


# Going to try removing the charsets block from the nexus file, and using the phylip file as the input alignment file, and using the nexus_for_SWSC.nexus file
# Turned it into a .txt file with just the charsets block, replaced charsts with "DNA" and its called 5_partitionfile_DNA.txt
# This qsub is 5_modelfinder.qsub

[inovick@scc1 modelfinder]$ pwd
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder

#!/bin/bash -l
#$ -P mullenl
#$ -N modelfinder_1 #Job name
#$ -m bea #Sends emails about the job
#$ -M inovick@bu.edu #What email do you want?
#$ -j y #Joins standard output and error into a single file
#$ -o modelfinder_1.qlog #Here is the logfile
#$ -l h_rt=24:00:00 #How long we run it for
#$ -pe omp 16 #How many cores we use

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "=========================================================="

module load module load iqtree/2.2.2.6

iqtree2 -s concat_matrix_75pis_loci_phylip.phylip -p 5_partitionfile_DNA.txt -m TESTMERGEONLY -mset JC69,TN93,TNe,K80,F81,HKY,SYM,GTR -rcluster 10 -nt AUTO

# This merged my partitions

# My partitions and recommended models are here
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/5_partitionfile_DNA.txt.best_scheme.nex

## SWSC-EN: good partitioning scheme for UCE data

# Script: 
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/PFinderUCE-SWSC-EN/py_script/SWSCEN.py

# File: 
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/concat_matrix_75pis_loci_nexus.nexus

# Command:
python /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/PFinderUCE-SWSC-EN/py_script/SWSCEN.py /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/nexus_for_SWSC.nexus

# Output: made an output in .cfg format called
nexus_for_SWSC.nexus_entropy_partition_finder.cfg


# Also found in 
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/nexus_for_SWSC.nexus_entropy_partition_finder.cfg

# Now going to try and merge partitions from the swsc .cfg output in modelfinder using the concat_matrix_75pis_loci_phylip.phylip and new dna partition text file, modified from nexus_for_SWSC.nexus_entropy_partition_finder.cfg, called SWSC_partitionfile_dna.txt
# Going to make a new qsub called swsc_partitionfinder_partitioning.qsub

iqtree2 -s concat_matrix_75pis_loci_phylip.phylip -p SWSC_partitionfile_dna.txt -m TESTMERGEONLY -mset JC69,TN93,TNe,K80,F81,HKY,SYM,GTR -rcluster 10 -nt AUTO

# My partitions and recommended models are here
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/SWSC_partitionfile_dna.txt.best_scheme.nex


# Now need to make nexus file for beauti on local machine by copying
concat_matrix_75pis_loci_phylip.phylip

# and then adding on the recommended models and partitions from 
SWSC_partitionfile_dna.txt.best_scheme.nex
# and 
5_partitionfile_DNA.txt.best_scheme.nex

# To make the final input files:
5_partitionfile_modelfinder_75pis_beast_supermatrix.nex #modelfinder partitioned data
swsc_partitionfile_modelfinder_inputforbeast_nexus.nex

## Mesquite:
# Making a new directory in beast called mesquite, making a new raxml tree with the 75 most informative sites to be used in mesquite
# Mesquite dir found here: 
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/mesquite

# Running raxml one more time on the subsetted 75 parsimonious informative loci on this file: 
concat_matrix_75pis_loci_phylip.phylip

# Found here: 
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/subsetting_loci/concat_matrix_75pis_loci_phylip

########## 1st BEAST run-- unlinked clocks birth-death 2/7/24
# Using BEAUTI1
# Following Ryan St. Laurent's instructions found in BEAST instructions file on local computer
# Using /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/5_partitionfile_modelfinder_70pis_beast_supermatrix.nex
# This partition file is made from modelfinder in iqtree2 using the following parameters (merged partitions based on modelfinders algorithm):

iqtree2 -s concat_matrix_75pis_loci_phylip.phylip -p 5_partitionfile_DNA.txt -m TESTMERGEONLY -mset JC69,TN93,TNe,K80,F81,HKY,SYM,GTR -rcluster 10 -nt AUTO

# And is a modified version of this file:

/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/5_partitionfile_DNA.txt.best_scheme.nex

# It uses the partitions and models from the above file, and uses the alignments from this file:

/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/modelfinder/concat_matrix_75pis_loci_nexus.nexus


# Partitions tab:
# Site models are unlinked, clock models are unlinked, trees are linked
# Rename tree partition model to default when clicking "link trees"

# Taxa tab:
# Taxon sets:
  # Cydia taxon set: contains cydia_amp and tineola_bisselliella_21_02
  # Adela taxon set: contains adela_rea and tineola_bisselliella_21_02

# Skip tips tab
# Skip traits tab

# Sites tab:
# Putting in nucleotide substitution models from supermatrix.nex file
# F= empirical state frequency so select empirical for those
# K2p= HKY
# If no F specified, it's equal base frequencies
# TN= TN93
# SYM=GTR with base frequencies set to all equal

# Clocks tab:
# Since unlinked, have to do a clock model for each partition. When doing linked, there is only one clock model
# I'm going to use the same clock model for each partition: uncorrelated relaxed clock model with lognormal relaxed distribution

# Trees tab:
# Tree prior: Speciation birth death process analyses, used for all subsets. When doing Yule, select Yule process
# Select random starting tree for all subsets

# Skip states tab
# Priors tab:
# Adela: make prior distribution uniform, upper=211.89, lower=172.18 (Kawahara et al 2019 secondary calibrations from supplementary tree)
# Cydia: make prior distribution uniform, upper=153.94, lower=122.61 (Kawahara et al 2019 secondary calibrations from supplementary tree)
# Tree model root height: prior distribution uniform, upper=229.37, lower=188.52 (Kawahara et al 2019 secondary calibrations from supplementary tree)
  
# Operators tab:
# Select "fixed tree topology" operator mix

# MCMC tab:
# Length of chain to 200 million
# Echo state to screen every 20,000
# Log paramteres every 20,000
# Path sampling stepping stone settings should be default settings

# Now also using fixed tree topology, done in mesquite 3.81 on old computer with java 8, used subsetted 75 loci raxml tree.
# xml saved here: /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/modelfinder_75pis_linkedclocks_birthdeath/v4_ryan_partitionfile_modelfinder_secondarycal_75pis_birthdeath_supermatrix.xml
# 2 xml files submitted: one is beast_modelfinder_75pis_linkedclocks_birthdeath_startingtree_secondarycals_fixedtree_rep1 that has opterator mix "fixed tree topology", and the other is beast_modelfinder_75pis_linkedclocks_birthdeath_startingtree_secondarycals_rep1, which has "classic operator mix". Both of these have a fixed tree in the xml


# For adding the fixed tree to the xml files, make sure you turn on operator mix "fixed tree topology" when creating xml files and then paste this into the xml file above line 5530 :
# Also make sure to check if this run is yule or birth death and modify accordingly
<!-- Using birth-death model on tree: Gernhard T (2008) J Theor Biol, Volume 253, Issue 4, Pages 769-778 In press-->
<birthDeathModel id="birthDeath" units="years">
<birthMinusDeathRate>
  <parameter id="birthDeath.meanGrowthRate" value="2.0" lower="0.0"/>
</birthMinusDeathRate>
<relativeDeathRate>
  <parameter id="birthDeath.relativeDeathRate" value="0.5" lower="0.0"/>
</relativeDeathRate>
</birthDeathModel>

<!-- This is a simple constant population size coalescent model              -->
<!-- that is used to generate an initial tree for the chain.                 -->
<constantSize id="initialDemo" units="substitutions">
    <populationSize>
        <parameter id="initialDemo.popSize" value="100.0"/>
    </populationSize>
</constantSize>
<newick id="startingTree">
((cydia_amp:150.0,(((opogona_nipponica_21_31:0.10861680221993808,opogona_nipponica_21_32:0.10861680221993808)100:25.272700547541337,(((erechthias_zebrina_22_26:1.6051251410569343,erechthias_zebrina_22_16:1.6051251410569343)100:22.123134978293578,((scardiella_approximatella:16.200595383928214,morophaga_bucephala:16.200595383928214)100:6.557785482034113,trichophaga_tapetzella:22.758380865962327)69:0.9698792533881857)61:0.6116843271460866,(((((monopis_pavlovski:0.3154038179252723,(monopis_monachella_21_26:0.2340438193683195,monopis_monachella_21_27:0.2340438193683195)68:0.08135999855695278)100:8.16259152416316,((monopis_marginistrigella:4.408408052790404,(mon_laev_genome:0.20968012239056755,monopis_laevigella:0.20968012239056755)100:4.198727930399836)100:2.090885592613825,(monopis_crocicapitella_21_22:0.2767464797392316,monopis_crocicapitella_22_04:0.2767464797392316)100:6.222547165664997)100:1.9787016966842037)100:8.116116752992513,(tin_semi_genome:10.544681983428587,(((phereoeca_spp_22_08:0.07554344361693532,phereoeca_spp_22_07:0.07554344361693532)100:1.234417885985525,(tinea_pellionella_22_01:0.17243655788025425,tinea_pellionella_22_02:0.17243655788025425)100:1.1375247717222061)100:7.928685326624006,((tinea_apicimaculella_21_12:5.686133072087592,(((tinea_pellionella_22_03:0.10636776168065604,niditinea_fuscella_22_10:0.10636776168065604)95:0.034200034178237304,niditinea_fuscella_22_09:0.14056779585889334)100:4.086223291799343,tin_trino_genome:4.226791087658237)100:1.4593419844293551)100:1.7055835497198863,(((tin_pell_genome:1.2892504767431436,tinea_murariella:1.2892504767431436)100:2.6629266679341246,tinea_occidentalis:3.952177144677268)100:0.5142760607414587,(tineola_bisselliella_21_02:0.41461941588780604,(tineola_bisselliella_22_24:0.28623458865549706,tineola_bisselliella_22_17:0.28623458865549706)100:0.12838482723230898)100:4.051833789530921)100:2.925263416388751)100:1.846930034418988)100:1.306035327202121)100:6.049430111652358)100:1.554913013584489,((doleromorpha_porphyria_22_20:0.15947188612067364,doleromorpha_porphyria_22_24:0.15947188612067364)100:4.687309105821743,tinea_carnariella:4.846780991942417)100:13.302244116723017)100:1.8046625393199278,((xystrologa_wielgusi_21_37:0.03682721124894783,xystrologa_wielgusi_21_38:0.03682721124894783)100:10.820113287146958,hybroma_servulella:10.856940498395906)100:9.096747149589456)81:4.386256798511237)74:1.0413729032646764)83:12.118682650238725,acrolophus_spp:37.5)83:112.5):50.0,adela_rea:200.0):35.0;
    </newick>
<treeModel id="treeModel">
    <newick idref="startingTree"/>
    <rootHeight>
        <parameter id="treeModel.rootHeight"/>
    </rootHeight>
    <nodeHeights internalNodes="true">
        <parameter id="treeModel.internalNodeHeights"/>
    </nodeHeights>
    <nodeHeights internalNodes="true" rootNode="true">
        <parameter id="treeModel.allInternalNodeHeights"/>
    </nodeHeights>
</treeModel>

##################Running beast

#Qsub headers and running prerequisites:

#!/bin/bash -l

#$ -P mullenl
#$ -N partitionfinder_1clock_birthdeath_rep1 #Job name
#$ -m bea #Sends emails about the job
#$ -M inovick@bu.edu #What email do you want?
#$ -j y #Joins standard output and error into a single file
#$ -o qlog #Here is the log directory
#$ -l h_rt=360:00:00 #How long we run it for
#$ -pe omp 16 #How many cores we use

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "=========================================================="

module purge
module use /projectnb/mullenl/software/module #Has to be in this structure: module/link to module file
module load java/17.0.8
module load beast/1.10.4
module load beagle-lib/4.0.1
module list


#Beast runs can be found here :
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/

######The modelfinder partitioned, 75pis loci, linked clocks, birth-death model with 3 replicates found in:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/modelfinder_75pis_linkedclocks_birthdeath
### Rep 1:
# qsub: modelfinder_75pis_linkedclocks_birthdeath.qsub modelfinder_75pis_linkedclocks_birthdeath_rep1.qsub
# xml: final_partitionfinder_1clock_birthdeath.xml

### Rep2:
# qsub: modelfinder_75pis_linkedclocks_birthdeath_rep2.qsub
# xml: final_partitionfinder_1clock_birthdeath.xml

### Rep3:
# qsub: modelfinder_75pis_linkedclocks_birthdeath_rep3.qsub
# xml: final_partitionfinder_1clock_birthdeath.xml

###### The modelfinder partitioned, 75pis loci, linked clocks, yule model with 3 replicates can be found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/modelfinder_75pis_linkedclocks_yule
### Rep1:
# qsub: modelfinder_75pis_linkedclocks_yule_rep1.qsub
# xml: final_partitionfinder_1clock_yule.xml

### Rep2:
# qsub: modelfinder_75pis_linkedclocks_yule_rep2.qsub
# xml: final_partitionfinder_1clock_yule.xml

### Rep3:
# qsub: modelfinder_75pis_linkedclocks_yule_rep3.qsub
# xml: final_partitionfinder_1clock_yule.xml

###### The modelfinder partitioned, 75pis loci, unlinked clocks, birth death model with 3 replicates can be found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/modelfinder_75pis_linkedclocks_birthdeath
### Rep1:
# qsub: modelfinder_75pis_unlinkedclocks_birthdeath_rep1.qsub
# xml: final_partitionfinder_unlinkedclocks_birthdeath.xml

### Rep2:
# qsub: modelfinder_75pis_unlinkedclocks_birthdeath_rep2.qsub
# xml: final_partitionfinder_unlinkedclocks_birthdeath.xml

### Rep3:
# qsub: modelfinder_75pis_unlinkedclocks_birthdeath_rep3.qsub
# xml: final_partitionfinder_unlinkedclocks_birthdeath.xml

###### The modelfinder partitioned, 75pis loci, unlinked clocks, yule model with 3 replicates can be found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/modelfinder_75pis_unlinkedclocks_yule
### Rep1:
# qsub: modelfinder_75pis_unlinkedclocks_yule_rep1.qsub
# xml: final_partitionfinder_unlinkedclocks_yule.xml

### Rep2:
# qsub: modelfinder_75pis_unlinkedclocks_yule_rep2.qsub
# xml: final_partitionfinder_unlinkedclocks_yule.xml

### Rep3:
# qsub: modelfinder_75pis_unlinkedclocks_yule_rep3.qsub
# xml:final_partitionfinder_unlinkedclocks_yule.xml

###### The SWSC partitionded, 75pis loci, linked clocks, birth death model with 3 replicates can be found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/swsc_75pis_linkedclocks_birthdeath
### Rep1:
# qsub: swsc_linkedclocks_birthdeath_rep1.qsub
# xml: final_swsc_1clock_birthdeath.xml

### Rep2:
# qsub: swsc_linkedclocks_birthdeath_rep2.qsub
# xml: final_swsc_1clock_birthdeath.xml

### Rep3:
# qsub: swsc_linkedclocks_birthdeath_rep3.qsub
# xml: final_swsc_1clock_birthdeath.xml

###### The SWSC partitionded, 75pis loci, linked clocks, yule model with 3 replicates can be found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/swsc_75pis_linkedclocks_yule
### Rep1:
# qsub: swsc_75pis_linkedclocks_yule_rep1.qsub
# xml: final_swsc_1clock_yule.xml

### Rep2:
# qsub: swsc_75pis_linkedclocks_yule_rep2.qsub
# xml: final_swsc_1clock_yule.xml

### Rep3:
# qsub: swsc_75pis_linkedclocks_yule_rep3.qsub
# xml: final_swsc_1clock_yule.xml

###### The SWSC partitionded, 75pis loci, unlinked clocks, birth death model with 3 replicates can be found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/swsc_75pis_unlinkedclocks_birthdeath
### Rep1:
# qsub: swsc_75pis_unlinkedclocks_birthdeath_rep1.qsub
# xml: final_swsc_unlinkedclocks_birthdeath.xml

### Rep2:
# qsub: swsc_75pis_unlinkedclocks_birthdeath_rep2.qsub
# xml: final_swsc_unlinkedclocks_birthdeath.xml

### Rep3:
# qsub: swsc_75pis_unlinkedclocks_birthdeath_rep3.qsub
# xml: final_swsc_unlinkedclocks_birthdeath.xml

###### The SWSC partitionded, 75pis loci, unlinked clocks, yule model with 3 replicates can be found here:
/projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/beast/beast/beastruns/swsc_75pis_unlinkedclocks_yule
### Rep1:
# qsub: swsc_75pis_unlinkedclocks_yule_rep1.qsub
# xml: final_swsc_unlinkedclocks_yule.xml

### Rep2:
# qsub: swsc_75pis_unlinkedclocks_yule_rep2.qsub
# xml: final_swsc_unlinkedclocks_yule.xml

### Rep3:
# qsub: swsc_75pis_unlinkedclocks_yule_rep3.qsub
# xml: final_swsc_unlinkedclocks_yule.xml



# After they have all run, you are interested in the .log, mle.result.log file and the .trees file
# Cyberduck to your computer

##### For log files:

##### logcombiner:
# Open logcombiner to combine log files-- all 3 reps from a single run type. Make sure they are all just ".log"
  # Add plus sign in bottom left corner, make sure file type selected is log File
  # Burnin: 50million (a quarter of all steps)
  # Name the output file "[that run] combined logs"

##### Tracer:
# Then open the combined log file in tracer using the plus sign under "trace files"
    # Want ESS values for joint, prior and likelihood to be above 200
    # Look at the trace (top left tab) to see the states. It should look big and zig zaggy, likelihoods should be up and down repeatedly

##### For tree files:
# logcombiner, select the 3 tree files for each replicate for each per run
# This is selecting 10,000 trees from each analysis. These combined trees then go in treeannotator
# Burnin: 2500, a quarter of 10000 trees
# This is really big and will take like 10 min
# Save as "[that run] combined trees"
# After combining tree files, delete the 3 original tree files from your computer because they are huge and will take up space


##### Tree annotator:
# Don't specify burnin
# Defaults for everything: max clade credibility tree etc.
# Input: choose combined tree file
# Output: "[That run] MCC tree"
# This also takes a long time to run
# Then look at it in figtree, select node labels: node ages

# Then look at the mle.result.log file for each rep for each run, and record the path sampling marginal likelihood and stepping stone sampling for each rep for each run


###############################################################################
#Summary of steps of UCE pipeline
1. Download the data (they are in fastq.gz format)
2. Count read data. This counts the actual number of reads in a given sequence file for each species. I had them in R1 and R2 files
#i7 is for R1 file, forward 5' to 3', i5 is for R2 file, is reverse complement of given i5 index
3. Clean the read data for adapter contam and low qual bases using illumiprocessor, giving it a config file that contains which adapters are in which R1 and R2 file. Indexes are 10 nucleotides long
#These got put into trimmed_fastq in trimming_dir
#Split-adapter-quality-trimmed directory has all the cleaned reads
4. Assemble the data: use spades. Input clean read data from Split-adapter-quality-trimmed. Output goes into /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/array_assembly/spades-assemblies-array/contigs
#These are species specific assembly files
#The contigs directory is the important one, because it contains symlinks to all of the species- specific contigs. This means that you can treat this single folder as if it contains all of your assembled contigs.
#Assembly QC: I have the outputs pasted above
5. Finding UCE loci: find the contigs which are UCE loci and filter out those that are not. Need the probe set used for enrichments (ours was the Lepidooptera set)
#Lepidoptera set has 14,363 baits targeting 1,381 conserved loci
#Run phyluce-assembly-match-contigs-to-probes program
#Log file with output found in /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/match_contigs_probes_dir/phyluce_assembly_match_contigs_to_probes.log
#This shows the capture data for each sample
#The search we just ran created lastz search result files for each taxon, and stored summary results of these searches in the probe.matches.sqlite database
6. Extract UCE loci. We need to choose which taxa we want in the analysis, create a list (cofig file) of the taxa, then make a "data matrix configuration file" that has a list of which UCE loci we enriched in each taxon. This is used to extract FASTA data for each taxon for each UCE locus
#Make the initial list of loci for each taxon using phylue-assembly-get-match-counts which searches the sqlite database for probe matches
#Log file for this found in /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/extract_uce_dir/2_taxon_sets/2_all/phyluce_assembly_get_match_counts.log. Output is 3-all-taxa-incomplete.conf
7. Extract FASTA data that correspond to the loci in 3-all-taxa-incomplete
#Uses phyluce-assembly-get-fastas-from-match-counts
#Output is 3-all-taxa-incomplete.fasta. Log is in the log directory and is called phyluce_assembly_get_fastas_from_match_counts.log
#Incomplete matrix is in 3-all-taxa-incomplete.incomplete
8. Exploding the monolithic FASTA file
#This gives the individual stats on UCE assemblies for a given taxon, and we can do this by exploding the jumbled fasta file into a file of UCE loci that we have enriched by taxon and then running stats on those separated files
#These are found in the exploded-fastas directory
9. Aligning UCE loci
#I edge trimmed and internal trimmed because my taxa presumably span a wide range of divergence times (>50MYA)
?#What does this mean ->"In phyluce, we implement our edge-trimming algorithm by running the alignment program “as-is” (i.e., without the –no-trim) option. We do internal-trimming by turning off trimming using –no-trim, then passing the resulting alignments (in FASTA format) to a parallel wrapper around Gblocks."
#I used the mafft aligner
10. Edge trimming
#Used phyluce-align-seqcap-align, input was 3-all-taxa-incomplete.fasta, output was 3-mafft-nexus-edge-trimmed. Log file is found in phyluce-align-seqcap-align.log
?#This spits out a list of nexus files for (I think) each UCE locus that was captured, and in the nexus file is a list of taxa that share that locus and their sequence alignments against each other?
#Summary stats of alignments log file found in log directory, called phyluce_align_get_align_summary_data.log
#Most important data is number of loci we have and number of loci in data matrices of different completeness

[Matrix 50%]            741 alignments
[Matrix 55%]            693 alignments
[Matrix 60%]            597 alignments
[Matrix 65%]            477 alignments
[Matrix 70%]            361 alignments
[Matrix 75%]            236 alignments
[Matrix 80%]            110 alignments
[Matrix 85%]            39 alignments
[Matrix 90%]            17 alignments
[Matrix 95%]            2 alignments

[Alignments] loci:	1,171

11. Internal trimming: running internal trimming on the resulting alignments
#input 3-all-taxa-incomplete.fasta
#output 3-mafft-nexus-internal-trimmed
?#spits out another list of fasta files instead of nexus, but in his tutorial he has them in nexus format even though he specific in the qsub that the output format should be in fasta. whaaaa?
12. Trim the loci using gblocks
#This uses phyluce-align-get-gblocks-trimmed-alignments-from-untrimmed
#alignments are in 3-mafft-nexus-internal-trimmed and the output is 3-mafft-nexus-internal-trimmed-gblocks
#This spits out a list of nexus files
13. Alignment cleaning
#The files in 3-mafft-nexus-internal-trimmed-gblocks are alignments that we have, and within each alignment file there are names that are combinations of taxa names and the locus name for that taxon, so we need to clean the alignmments because this is not what we want downstream apparently
#Use gblocks trimmed alignments, continues to use gblocks to clean Alignments
#Using phyluce-align-remove-locus-name-from-files
#input 3-mafft-nexus-internal-trimmed-gblocks and outputs 3-mafft-nexus-internal-trimmed-gblocks-clean
#I must have skipped this step because I don't actually have the above output file lol
14. Final data matrices
#I did a 75% completeness matrix-- means that in a study of 100 taxa total, all alignments will contain at least 75 of these 100 taxa
#output is 3-mafft-nexus-internal-trimmed-gblocks-clean-75p
#I have 236 alignments of 1171 loci at 75% completeness
#Inside that directory are nexus files
#Looking at the log file phyluce-align-get-only-loci-with-min-taxa.log

Copied 236 alignments of 1171 total containing ≥ 0.75 proportion of taxa (n = 28)

15. Prepping for downstream analysis
#This is when we turn the 75p data that are in nexus files into a phylip File
#What is in the each nexus file are individual locus alignments for each taxon at that locus.
#What is in the final phylip or nexus file is a concatenation of the alignments


#Use phyloch for parsimony informative sites, after running script, you will have file called inform_names_PIS_15.txt
#Go into 3-mafft-nexus-internal-trimmed-gblocks-clean-75p and copy the loci from the list into a new directory
#Then concat into single alignment
#Then run beast using that


##### For publication
# Data found here: /projectnb/mullenl/novick/uce_analysis/uce-analysis-pipeline/trimming_dir/uce-data
