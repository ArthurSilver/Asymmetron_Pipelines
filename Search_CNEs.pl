####align genome using lastz##############
#make_lastz_chains/make_chains.py Al Hs Al.rename.fasta GRCh38.rename.fasta --pd Al_Hs -f --chaining_memory 8
#chainToAxt Al.Hs.final.chain Al.rename.2bit Hs.rename.2bit Al.Hs.final.axt

####record the chromosome length infomation
open(DATA,"</$work_dir/CNE/GRCh38.rename.fasta.fai") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom_len{$temp[0]}{"Hs"}=$temp[1];
}
close(DATA);

open(DATA,"</$work_dir/CNE/Mm.rename.fasta.fai") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom_len{$temp[0]}{"Mm"}=$temp[1];
}
close(DATA);

open(DATA,"</$work_dir/CNE/Oa.rename.fasta.fai") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom_len{$temp[0]}{"Oa"}=$temp[1];
}
close(DATA);

open(DATA,"</$work_dir/CNE/Lo.rename.fasta.fai") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom_len{$temp[0]}{"Lo"}=$temp[1];
}
close(DATA);

open(DATA,"</$work_dir/CNE/Xt.rename.fasta.fai") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom_len{$temp[0]}{"Xt"}=$temp[1];
}
close(DATA);

open(DATA,"</$work_dir/CNE/Gg.rename.fasta.fai") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom_len{$temp[0]}{"Gg"}=$temp[1];
}
close(DATA);


####extract he block which length >20 and identify>60%, and convert the coordinates of minus strand
@file=glob("/$work_dir/CNE/axt/*.axt");

for $filename(@file){
	$filename=~/axt\/(.*)\.axt/;
	my %seq;
	my %record;
	print "processing $1\n";
	open(INPUT,">$1.bed");
	open(DATA,"<$filename");
	while(<DATA>){
		chomp;
		$filename=~/axt\/Alref\.(.*)\.final.axt/;
		$sp=$1;
		if(/^\d/){
			@temp=split(/\s+/);
			$key="$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]";
			if($temp[7] eq "-"){
				die "can't match chrom length" if(!exists($chrom_len{$temp[4]}{$sp}));
				$start=$chrom_len{$temp[4]}{$sp}-$temp[6]+1;
				$end=$chrom_len{$temp[4]}{$sp}-$temp[5]+1;
				$key="$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$start\t$end\t$temp[7]";
			}
			$record{$key}=0;
			$count=1;
			next;
		}
		else{
			$seq{$key}{$count}="$_";
			$count++;
		}	
	}
	
	for $key(sort keys %record){
		@seq1=split(//,$seq{$key}{1});
		@seq2=split(//,$seq{$key}{2});
		$len=length($seq{$key}{1});
		my $count=0;
		for $i(0..$#seq1){
			next if($seq1[$i] eq "-" || $seq2[$i] eq "-");
			if(uc($seq1[$i]) eq uc($seq2[$i])){$count++;}
		}
		$ratio=$count / $len;
		if($ratio>=0.6 && $len>=20){
			print INPUT "$key\t";
			printf INPUT "%.2f", $ratio;
			print INPUT "\t$len\n";
		}
	}
	
	close(INPUT);
	close(DATA);
}


######sort and merge bed file
@file=glob("/$work_dir/CNE/axt/Al*.axt");

for $filename(@file){
	$filename=~/axt\/(.*)\.axt/;
	$prefix=$1;
	print "processing $prefix\n";
	open(INPUT,">$prefix.tmp.bed");
	open(DATA,"<$prefix.bed");
	while(<DATA>){
		chomp;
		@temp=split(/\s+/);
		print INPUT "$temp[0]\t$temp[1]\t$temp[2]\n";
	}
	close(INPUT);
	close(DATA);
	system("sort -k1,1 -k2,2n $prefix.tmp.bed >$prefix.sort.bed");
	system("rm $prefix.tmp.bed");
	system("bedtools merge -i $prefix.sort.bed >$prefix.sort.merge.bed");
}


#####intersect target candidates using bedtools
#bedtools intersect -a Alref.Hs.final.sort.merge.bed -b Alref.Mm.final.sort.merge.bed >Alref.intersect2.bed
#
#bedtools intersect -a Alref.intersect2.bed -b Alref.Oa.final.sort.merge.bed >Alref.intersect3.bed
#
#bedtools intersect -a Alref.intersect3.bed -b Alref.Gg.final.sort.merge.bed >Alref.intersect4.bed
#
#bedtools intersect -a Alref.intersect4.bed -b Alref.Lo.final.sort.merge.bed >Alref.intersect5.bed
#
#bedtools intersect -a Alref.intersect5.bed -b Alref.Xt.final.sort.merge.bed >Alref.intersect6.bed
#
#bedtools intersect -a Alref.intersect6.bed -b Dp5.anno.block.len20.ide60.txt >Alref.intersect7.bed


#########convert Al coordinates to Hs and make pair###########

open(DATA,"</$work_dir/CNE/Hs.rename.fasta.fai") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom_len{$temp[0]}=$temp[1];
}
close(DATA);

my %seq;
my %record;
my %position;

open(DATA,"</$work_dir/CNE/axt/Alref.Hs.final.axt") or die "can't open data";
open(INPUT,">/$work_dir/CNE/axt/convert_coordinate.Hs.txt");
while(<DATA>){
	chomp;
	if(/^\d/){
		@temp=split(/\s+/);
		$key="$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]";
		if($temp[7] eq "-"){
			$start=$chrom_len{$temp[4]}-$temp[6]+1;
			$end=$chrom_len{$temp[4]}-$temp[5]+1;
			$key="$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$start\t$end\t$temp[7]";
		}
		$record{$key}=0;
		$count=1;
		next;
	}
	else{
		$seq{$key}{$count}="$_";
		$count++;
	}	
}
close(DATA);

my %filter;
open(DATA,"</$work_dir/CNE/axt/Alref.Hs.final.bed") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$filter{"$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]"}++;
}
close(DATA);

for $key(sort keys %record){
	next unless(exists($filter{$key}));
	@temp=split(/\s+/,$key);
	@seq1=split(//,$seq{$key}{1});
	@seq2=split(//,$seq{$key}{2});
	
	$a=0;$b=0;
	for $i(0..$#seq1){
		if($seq1[$i] ne "-"){$position_al=$temp[1]+$a;$a++;}
		if($temp[6] eq "+"){
			if($seq2[$i] ne "-"){$position_hs=$temp[4]+$b;$b++;}
		}
		elsif($temp[6] eq "-"){
			if($seq2[$i] ne "-"){$position_hs=$temp[5]-$b;$b++;}
		}
		if($seq1[$i] ne "-" && $seq2[$i] ne "-"){
			push(@{$position{$temp[0]}{$position_al}},"$temp[3]\t$position_hs");
		}
	}
}


open(DATA,"</$work_dir/CNE/axt/Alref.intersect7.bed") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	for $i($temp[1]..$temp[2]){
		for $position_hs (@{$position{$temp[0]}{$i}}){
			print INPUT "$temp[0]\t$i\t$position_hs\n";
		}
	}
}

my %coord_record;
open(DATA,"</$work_dir/CNE/axt/convert_coordinate.Hs.txt") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	push(@{$coord_record{$temp[0]}{$temp[1]}},"$temp[2]\t$temp[3]");
}
close(DATA);

open(DATA,"</$work_dir/CNE/axt/Alref.intersect7.bed");
open(INPUT,">/$work_dir/CNE/axt/candidate_pair.txt");

while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$len=$temp[2]-$temp[1]+1;
	next unless($len>=20);
	my %record;
	my %position;
	print "processing $_\n";
	for $i($temp[1]..$temp[2]){
		for $line(@{$coord_record{$temp[0]}{$i}}){
			my ($Hs_chrom,$Hs_coord)=split(/\s+/,$line);
			push(@{$record{$Hs_chrom}},$Hs_coord);
			$position{$Hs_chrom}{$Hs_coord}=$i;
		}
	}
	
	for $key(sort keys %record){
		my %hash=();
		$target_chrom=$key;
		
		@array=grep{++$hash{$_}<2} sort{$a<=>$b} @{$record{$key}};
		my $record_start="Null";
		my $record_position;
		for $i(0..$#array){
			if($record_start eq "Null"){
				if($array[$i+1]-$array[$i]<=5){
					$record_start=$i;
					$record_position=$i;
					next;
				}
				else{next;}
			}
			
			if($array[$i]-$array[$record_position]<=5){
				$record_position=$i;
			}
			
			else{
				print INPUT "$_\t";
				print INPUT "${target_chrom}\t$array[$record_start]\t$array[$record_position]\t";
				print INPUT "$temp[0]\t$position{$target_chrom}{$array[$record_start]}\t$position{$target_chrom}{$array[$record_position]}\n";
				
				if($array[$i+1]-$array[$i]<=5){
					$record_start=$i;
					$record_position=$i;
					next;
				}
				
				else{
					$record_start="Null";
				}
			}
			
			if($i==$#array && $record_start ne "Null"){
				if($record_position==$i){
					print INPUT "$_\t";
					print INPUT "${target_chrom}\t$array[$record_start]\t$array[$record_position]\t";
					print INPUT "$temp[0]\t$position{$target_chrom}{$array[$record_start]}\t$position{$target_chrom}{$array[$record_position]}\n";
				}
			}
		}
	}
	close(DATA2);
}
####################


open(DATA,"</$work_dir/CNE/axt/candidate_pair.txt");
open(INPUT,">/$work_dir/CNE/axt/Al_candidate.txt");
open(INPUT2,">/$work_dir/CNE/axt/Hs_candidate.txt");

while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$len1=$temp[5]-$temp[4]+1;
	if($temp[7]>$temp[8]){
		$start=$temp[8];
		$end=$temp[7];
		$ori="-";
	}
	else{
		$start=$temp[7];
		$end=$temp[8];
		$ori="+";
	}
	$len2=$end-$start+1;
	$diff=$len2-$len1;
	$diff=abs($diff);
	
	next if($len1<20 || $len2<20); #########sequence length >20 in two species
	if($diff<=20){
		print INPUT "$temp[3]\t";
		print INPUT $temp[4]-1;
		print INPUT "\t$temp[5]\n";
		print INPUT2 "$temp[6]\t";
		print INPUT2 $start-1;
		print INPUT2 "\t$end\t$ori\n";
	}
}
close(DATA);

############################


######################extract sequence from block pair, align again using muscle and the second screening###########
my @Hs_seq;my @Al_seq;
my %Hs_seq;my %Al_seq;
open(DATA,"<Hs_candidate.fasta");#####need to extract fasta sequence using bedtools
while(<DATA>){
	chomp;
	if(/>/){
		$key=$_;
		push(@Hs_seq,$key);
	}
	else{
		$Hs_seq{$key}=$_;
	}
}
close(DATA);

open(DATA,"<Al_candidate.txt");
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	push(@ori,$temp[-1]);
}

open(DATA,"<Al_candidate.fasta");
while(<DATA>){
	chomp;
	if(/>/){
		$key=$_;
		push(@Al_seq,$key);
	}
	else{
		if($ori[$#Al_seq] eq "+"){
			$Al_seq{$key}=$_;
		}
		else{
			@seq=split(//);
			@seq_reverse=reverse(@seq);
			my @seq_new;
			for $base(@seq_reverse){
				if($base eq "A" || $base eq "a"){push(@seq_new,"T");}
				if($base eq "T" || $base eq "t"){push(@seq_new,"A");}
				if($base eq "C" || $base eq "c"){push(@seq_new,"G");}
				if($base eq "G" || $base eq "g"){push(@seq_new,"C");}
			}
			$line=join("",@seq_new);
			$Al_seq{$key}=$line;
		}
	}
}
close(DATA);

##########output the pair sequence fasta file########
for $i(0..$#Hs_seq){
	open(INPUT,">/$work_dir/CNE/axt/seq/seq_pair/$i.fasta");
	print INPUT "$Hs_seq[$i]\n$Hs_seq{$Hs_seq[$i]}\n";
	print INPUT "$Al_seq[$i]\n$Al_seq{$Al_seq[$i]}\n";
}

for $i(0..$#Hs_seq){
	system("muscle -in seq_pair/$i.fasta -out seq_aln/$i.aln");
}

##############calculate length and indentity of pair alignment again#######
open(INPUT,">seq_pair_ratio.txt");
my %record;
for $i(0..$#Hs_seq){
	my %seq=();
	open(DATA,"</$work_dir/CNE/axt/seq/seq_aln/$i.aln") or die "can't open data";
	while(<DATA>){
		chomp;
		if(/>(.*)/){
			$key=$1;
		}
		else{
			$seq{$key}.=$_;
		}
	}
	close(DATA);
	$Hs_seq[$i]=~/>(.*)/;
	$key1=$1;
	@seq1=split(//,$seq{$1});
	
	$Al_seq[$i]=~/>(.*)/;
	$key2=$1;
	@seq2=split(//,$seq{$1});
	my $count=0;
	my $len=0;
	for $j(0..$#seq1){
		$len++;
		if($seq1[$j] eq $seq2[$j]){
			$count++;
		}	
	}
	$ratio=$count / $len;
	$record{$i}=$ratio;
	print INPUT "$key1\t$key2\t$len\t";
	printf INPUT "%.2f",$ratio;
	print INPUT "\n";
	print INPUT "$seq{$key1}\n$seq{$key2}\n"
}

##########screen the block pair which identity>60############
open(DATA,"</$work_dir/CNE/axt/seq_pair_ratio.txt");
open(INPUT,">/$work_dir/CNE/axt/Al_candidate_filter.txt");
open(INPUT2,">/$work_dir/CNE/axt/Hs_candidate_filter.txt");

while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	next unless($temp[-1]>=0.6);
	$temp[1]=~/(.*):(\d+)-(\d+)/;
	print INPUT "$1\t$2\t$3\n";
	$temp[0]=~/(.*):(\d+)-(\d+)/;
	print INPUT2 "$1\t$2\t$3\n";
}
close(DATA);
close(INPUT);


##############annotate the block pair and screen the results based on the annotation############

open(DATA,"</$work_dir/CNE/axt/tmp/Hs.gene.sort.gff") or die "can't open data"; #####annotation file of Human
my %record;
my %cds;
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom=$temp[0];
	$start=$temp[1];
	$end=$temp[2];
	next if($temp[3] eq "gene");
	if($temp[3] eq "CDS"){
		$cds{$chrom}{$start}{$end}="$temp[4]_$temp[3]";
		next;
	}
	$record{$chrom}{$start}{$end}="$temp[4]_$temp[3]";
	if($#temp==3){$record{$chrom}{$start}{$end}="$temp[3]";}
}
close(DATA);


open(DATA,"</$work_dir/CNE/axt/seq/Al_Hs/Hs_candidate.txt") or die "can't open data";
open(INPUT,">/$work_dir/CNE/axt/seq/Al_Hs/Hs_candidate.anno.txt");

while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$chrom=$temp[0];
	$start=$temp[1]+1;
	$end=$temp[2];
	my @record;
	print INPUT "$chrom\t$start\t$end\t";
	for $key(sort{$a<=>$b} keys %{$record{$chrom}}){
		for $key2(sort{$a<=>$b} keys %{$record{$chrom}{$key}}){
			if($start<=$key && $key2<=$end){
				push(@record,$record{$chrom}{$key}{$key2});
			}
			elsif($key<=$start && $start<=$key2){
				push(@record,$record{$chrom}{$key}{$key2});
			}
			elsif($key<=$end && $end<=$key2){
				push(@record,$record{$chrom}{$key}{$key2});
			}
		}
	}
	
	for $key(sort{$a<=>$b} keys %{$cds{$chrom}}){
		for $key2(sort{$a<=>$b} keys %{$cds{$chrom}{$key}}){
			if($start<=$key && $key2<=$end){
				push(@record,$cds{$chrom}{$key}{$key2});
			}
			elsif($key<=$start && $start<=$key2){
				push(@record,$cds{$chrom}{$key}{$key2});
			}
			elsif($key<=$end && $end<=$key2){
				push(@record,$cds{$chrom}{$key}{$key2});
			}
		}
	}
	
	my %hash;
	@record=sort grep{++$hash{$_}<2} @record;
	$line=join("\t",@record);
	print INPUT "$line\n";
}
close(DATA);
close(INPUT);


#########annotate the non genetic region, find their nearest gene###########
my %record;
sub min{
	@array=@_;
	$min=$array[0];
	for $i(1..$#array){
		if($array[$i]<$min){
			$min=$array[$i];
		}
	}
	return $min;
}
open(DATA,"</$work_dir/CNE/axt/tmp/Hs.gene.sort.gff") or die "can't open data";
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$record{$temp[0]}{"$temp[1]_$temp[2]"}="$temp[4]_$temp[3]";
}
close(DATA);

my %distance;
open(DATA,"</$work_dir/CNE/axt/seq/Al_Hs/Hs_candidate.anno.txt");
open(INPUT,">/$work_dir/CNE/axt/seq/Al_Hs/Hs_candidate.anno2.txt");
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	if($#temp>2){
		print INPUT "$_\n";
		next;
	}
	
	my %distance;
	for $position(sort keys %{$record{$temp[0]}}){
		my($start,$end)=split("_",$position);
		$dis1=abs($start-$temp[1]);
		$dis2=abs($end-$temp[1]);
		$dis3=abs($start-$temp[2]);
		$dis4=abs($end-$temp[2]);
		$dis=min($dis1,$dis2,$dis3,$dis4);
		$distance{$position}=$dis;
	}
	for $position(sort{$distance{$a}<=>$distance{$b}} keys %distance){
		print INPUT "$_\t";
		my($start,$end)=split("_",$position);
		if($distance{$position}<=100000){
			print INPUT "Near_$record{$temp[0]}{$position}\t$distance{$position}\n";
		}
		else{
			print INPUT "NearestGene>100k\n";
		}
		last;
	}
}

############annotate the Al candidates using same scripts###########

#######screen the block pair based on annotation (same genetic region and homologous genes)#############

open(DATA,"</$work_dir/CNE/Orthogroups.replace.tsv");
while(<DATA>){
	chomp;
	@temp=split(/\t/);
	next if($temp[1]!~/\w/ || $temp[7]!~/\w/);
	@Al_gene=split(/,/,$temp[1]);
	@Hs_gene=split(/,/,$temp[7]);
	for $al_gene(@Al_gene){
		for $hs_gene(@Hs_gene){
			$record{$al_gene}{$hs_gene}=$temp[0];
		}
	}
}
close(DATA);

open(DATA,"</$work_dir/CNE/Al_candidate.anno.txt");
open(INPUT,">/$work_dir/CNE/final_results.txt");
my $count=0;
my %position;my %anno;my %gene;
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$position{$count}{"Al"}="$temp[0]\t$temp[1]\t$temp[2]";
	
	if(/CDS/ && $_!~/Near/){
		$anno{$count}{"Al"}="CDS";
	}
	elsif(/NearestGene>100k/){
		$anno{$count}{"Al"}="noAnno";
	}
	else{
		if(/Near/){
			$anno{$count}{"Al"}="Intergenic";
		}
#		elsif(/intron/){
#			$anno{$count}{"Al"}="Intron";
#		}
		else{
			$anno{$count}{"Al"}="NonCoding";
		}
	}
	
	for $i(3..$#temp){
		if($temp[$i]=~/(Alref_\d+)_/){
			$gene{$count}{"Al"}{$1}++;
		}
	}
	$count++;
}
close(DATA);

open(DATA,"</$work_dir/CNE/Hs_candidate.anno.txt");
my $count=0;
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$position{$count}{"Hs"}="$temp[0]\t$temp[1]\t$temp[2]";
	
	if(/CDS/ && $_!~/Near/){
		$anno{$count}{"Hs"}="CDS";
	}
	elsif(/NearestGene>100k/){
		$anno{$count}{"Hs"}="noAnno";
	}
	
	else{
		if(/Near/){
			$anno{$count}{"Hs"}="Intergenic";
		}
#		elsif(/intron/){
#			$anno{$count}{"Hs"}="Intron";
#		}
		else{
			$anno{$count}{"Hs"}="NonCoding";
		}
	}
	
	for $i(3..$#temp){
		if($temp[$i]=~/(gene-.*?)_/){
			$gene{$count}{"Hs"}{$1}++;
		}
	}
	$count++;
}
close(DATA);

for $i(0..$count-1){
	$record_gene="";
	next unless($anno{$i}{"Hs"} eq $anno{$i}{"Al"});
	next if($anno{$i}{"Hs"} eq "CDS");
	if($anno{$i}{"Hs"} eq "noAnno"){
		print INPUT "$position{$i}{'Al'}\t$position{$i}{'Hs'}\t$anno{$i}{'Hs'}\n";
		next;
	}
	my $gene_count=0;
	for $Al_gene(sort keys %{$gene{$i}{"Al"}}){
		for $Hs_gene(sort keys %{$gene{$i}{"Hs"}}){
			if(exists($record{$Al_gene}{$Hs_gene})){
				$gene_count++;
				$record_gene="$Al_gene\t$Hs_gene";
			}
		}
	}
	next unless($gene_count);
	print INPUT "$position{$i}{'Al'}\t$position{$i}{'Hs'}\t";
	print INPUT "$anno{$i}{'Al'}\t$anno{$i}{'Hs'}\t$record_gene\n";	
}