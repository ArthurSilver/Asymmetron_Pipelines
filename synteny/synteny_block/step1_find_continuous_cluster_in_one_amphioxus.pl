#!/usr/bin/perl
##############################################################################
# Script Name: step1_find_continuous_cluster_in_one_amphioxus.pl
# Purpose: Identify continuous gene clusters in amphioxus genome 
#          using gene position data from MCscanX results
# Input Files: 
#   1. Alref.gff - Contains gene locus information (chromosome, start, end, gene ID)
#   2. 1.txt     - Sorted list of gene IDs (sorted by chromosome and position)
# Output Files:
#   1. 2.txt     - Gene IDs with positional distance to previous gene
#   2. 3.txt     - Gene IDs annotated with cluster ID and orientation/start status
#   3. 4.txt     - Final output with filtered clusters (only clusters ¡Ý3 genes retained)
# Notes: This script must be run separately for each of the five amphioxus species
##############################################################################

# --------------------------
# Step 1: Parse GFF file to map genes to chromosomes and store gene order
# --------------------------
# Open GFF file containing gene locus information (format: Chr1 27 3225 Alref_000001-T1)
open(DATA,"<file_path/Alref.gff");  

# %chrom: Hash to store chromosome ID for each gene (key=gene ID, value=chromosome)
my %chrom;
# %record: Hash to store ordered list of genes per chromosome (key=chromosome, value=array of gene IDs)
my %record;

# Read GFF file line by line
while(<DATA>){
    chomp;  # Remove newline character
    @temp = split(/\s+/);  # Split line by whitespace into array
    
    # Map gene ID (4th field) to its chromosome (1st field)
    $chrom{$temp[3]} = $temp[0];
    # Add gene ID to the chromosome's gene list in %record
    push(@{$record{$temp[0]}}, $temp[3]);
}
close(DATA);  # Close GFF file

# --------------------------
# Step 2: Calculate positional distance between consecutive genes
# --------------------------
my $count = 0;  # Counter to track first line in input file
# Open sorted gene ID file (1.txt) and output file (2.txt) for distance values
open(DATA,"<file_path/1.txt");  
open(INPUT,">file_path/2.txt");

# Read sorted gene IDs line by line
while(<DATA>){
    chomp;
    if($count == 0){  # Handle first line (no previous gene to compare)
        print INPUT "$_\t0\n";  # Assign distance 0 to first gene
        $record = $_;  # Store first gene ID as previous gene
        $count++;
        next;
    }
    else{  # Process subsequent genes
        my ($a, $b) = ($record, $_);  # $a=previous gene, $b=current gene
        my $chr_a = $chrom{$a};       # Chromosome of previous gene
        my $chr_b = $chrom{$b};       # Chromosome of current gene
        
        # Get full gene list for chromosomes of previous/current gene
        my @array_a = @{$record{$chr_a}};
        my @array_b = @{$record{$chr_b}};
        
        # Find index position of $a in its chromosome's gene list
        my @tmp_a = grep { $array_a[$_] eq $a } 0 .. $#array_a;
        # Find index position of $b in its chromosome's gene list
        my @tmp_b = grep { $array_b[$_] eq $b } 0 .. $#array_b;
        
        if($chr_a ne $chr_b){  # Genes on different chromosomes
            print INPUT "$_\t100\n";  # Assign large distance (100) for cross-chromosome
        }
        else{  # Genes on same chromosome
            $dis = $tmp_b[0] - $tmp_a[0];  # Calculate index distance between genes
            print INPUT "$_\t$dis\n";      # Write current gene + distance
        }
        
        $record = $_;  # Update previous gene to current gene for next iteration
    }
}
close(DATA);
close(INPUT);

# --------------------------
# Step 3: Assign cluster IDs and orientation (plus/minus/start)
# --------------------------
open(DATA,"<file_path/2.txt");       # Open distance file (2.txt)
open(INPUT,">file_path/3.txt");      # Open output file for cluster annotations

my $count = 1;          # Cluster ID counter (starts at 1)
my $cluster_size = 0;   # Track number of genes in current cluster
my $ori;                # Store cluster orientation (plus/minus)

while(<DATA>){
    chomp;
    @temp = split(/\s+/);  # Split line into [gene ID, distance]
    
    if($temp[1] == 0){  # First gene (distance = 0)
        print INPUT "$_\tcluster$count\tstart\n";  # Mark as cluster start
        $cluster_size++;  # Increment cluster size
    }
    # Case 1: Positive distance (0 < distance ¡Ü5) - plus orientation
    elsif($temp[1] >= 0 && $temp[1] <=5){
        if($cluster_size == 1){  # First gene in new cluster
            $ori = "plus";
            print INPUT "$_\tcluster$count\t$ori\n";
            $cluster_size++;
        }
        else{  # Existing cluster
            if($ori eq "plus"){  # Consistent orientation - continue cluster
                print INPUT "$_\tcluster$count\t$ori\n";
                $cluster_size++;
            }
            else{  # Orientation change - start new cluster
                $cluster_size = 1;
                $count++;
                print INPUT "$_\tcluster$count\tstart\n";
            }
        }
    }
    # Case 2: Negative distance (-5 ¡Ü distance <0) - minus orientation
    elsif($temp[1] >= -5 && $temp[1] <0){
        if($cluster_size == 1){  # First gene in new cluster
            $ori = "minus";
            print INPUT "$_\tcluster$count\t$ori\n";
            $cluster_size++;
        }
        else{  # Existing cluster
            if($ori eq "minus"){  # Consistent orientation - continue cluster
                print INPUT "$_\tcluster$count\t$ori\n";
                $cluster_size++;
            }
            else{  # Orientation change - start new cluster
                $cluster_size = 1;
                $count++;
                print INPUT "$_\tcluster$count\tstart\n";
            }
        }
    }
    # Case 3: Distance outside ¡À5 - start new cluster
    elsif($temp[1] < -5 || $temp[1] >5){
        $cluster_size = 1;
        $count++;
        print INPUT "$_\tcluster$count\tstart\n";
    }
}
close(DATA);
close(INPUT);

# --------------------------
# Step 4: Filter clusters (retain only clusters with ¡Ý3 genes)
# --------------------------
open(DATA,"<file_path/3.txt");       # Open cluster-annotated file
open(INPUT,">file_path/4.txt");      # Final output file

# First pass: Count number of genes per cluster
my %count;
while(<DATA>){
    chomp;
    @temp = split(/\s+/);  # Split line into [gene ID, distance, cluster ID, status]
    $count{$temp[2]}++;    # Increment count for this cluster
}
close(DATA);

# Second pass: Filter clusters (mark small clusters as "No")
open(DATA,"<file_path/3.txt");
while(<DATA>){
    chomp;
    @temp = split(/\s+/);
    # If cluster has <3 genes, mark as "No" (non-cluster)
    if($count{$temp[2]} <3){
        $temp[2] = "No";
    }
    # Reconstruct line with tab separation
    $line = join("\t",@temp);
    print INPUT "$line\n";
}
close(DATA);
close(INPUT);

##############################################################################
# End of Script Notes:
# - This script processes one amphioxus species at a time
# - Requires separate execution for all five amphioxus species
# - Final output (4.txt) contains:
#   Column 1: Gene ID
#   Column 2: Distance to previous gene
#   Column 3: Cluster ID (or "No" for non-clusters)
#   Column 4: Orientation/start status (plus/minus/start)
##############################################################################