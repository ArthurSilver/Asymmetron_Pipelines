#######step1 split the maf file#########
use Getopt::Long;

GetOptions('target_chrom|t:s' => \$target_chrom,
           'work_dir|d:s' => \$work_dir); 

system("mkdir $work_dir/$target_chrom/phylop");
system("mkdir $work_dir/$target_chrom/stats");
system("mkdir $work_dir/$target_chrom/maf");

open(DATA,"<$work_dir/$target_chrom/Hs.$target_chrom.maf") or die "can't open data";
my %block;

my $count=0;
while(<DATA>){
	chomp;
	next unless($_);
	next if(/^#/);
	@temp=split(/\s+/);
	if(/^a/){
		$count=1;
		next;
	}
	@temp=split(/\s+/);
	if($count){
		$key1=$temp[1];
		$key2=$temp[2];
		$key3=$temp[3];
		$count=0;
	}
	push(@{$block{$key1}{$key2}{$key3}},$_);
}

my $len=0;
my $count=1;
open(INPUT,">$work_dir/$target_chrom/maf/Hs.$target_chrom.block1.maf");
for $key1(sort keys %block){
	for $key2(sort{$a<=>$b} keys %{$block{$key1}}){
		for $key3(sort keys %{$block{$key1}{$key2}}){
#			next unless($key3>=10);
			next if($#{$block{$key1}{$key2}{$key3}}==0);
			my @species_count;
			for $line(@{$block{$key1}{$key2}{$key3}}){
				@temp=split(/\s+/,$line);
			  $species=(split(/\./,$temp[1]))[0];
			  if($species eq "Alref" || $species eq  "Bb" || $species eq "Bf" || $species eq "Bj" || $species eq "Bl" || $species eq "Gg" || $species eq "Hs" || $species eq "Xt" || $species eq "Lo" || $species eq "Oa" || $species eq "Mm"){
				  push(@species_count,$species);
			  }
			}
			my %hash;
			@species_count=grep{++$hash{$_}<2} @species_count;
			next if($#species_count==0);
			next unless($#species_count==10);
			$len_record=$len;
			$len+=$key3;
			if($len_record>=1000000){
				print INPUT "##block total length is $len_record\n";
				$len=0;
				$count++;
				close(INPUT);
				open(INPUT,">$work_dir/$target_chrom/maf/$key1.block$count.maf"); 
			}
			print INPUT "a\n";
			my $Hs_count=0;
			for $line(@{$block{$key1}{$key2}{$key3}}){
				@temp=split(/\s+/,$line);
				$species=(split(/\./,$temp[1]))[0];
				if($species eq "Hs"){
					$Hs_count++;
					next if($Hs_count>=2);					
				}
				next if($species=~/A\d+/);
				print INPUT "$line\n";
			}
		}
	}
}

##########step2 calculate phylop############
@file=glob("$work_dir/$target_chrom/maf/*.maf");

for $filename(@file){
	$filename=~/(Hs.NC.+\.block\d+).maf/;
	$id=$1;
	print "processing $id\n";
	system("phyloP --msa-format MAF --wig-scores --method LRT --mode CONACC $work_dir/chordata.mod $filename >$work_dir/$target_chrom/phylop/$id.wig");
}

###########step3 count the depth and phylop of bases#####################
@file=glob("$work_dir/$target_chrom/maf/*.maf") or die "can't open glob";

for $filename(@file){
	$filename=~/Hs.(NC.*?)\.(block\d+).maf/;
	$chrom="$1";$block=$2;
	open(DATA,"<$work_dir/$chrom/phylop/Hs.$chrom.$block.wig") or die "can't open phylop";
	my %phylop;
	while(<DATA>){
		chomp;
		if(/start=(\d+)/){
			$position=$1;
			next;
		}
		$phylop{$chrom}{$position}="$_";
		$position++;
	}
	
	close(DATA);
	
	open(DATA,"<$filename");
	my %len;
	my %seq;
	while(<DATA>){
		chomp;
		next if(/^#/);
		if(/^a/){
			$count=0;
			next;
		}
		if($count==0){
			@temp=split(/\s+/);
			$count++;
			$start=$temp[2]+1;
			$len{$chrom}{$start}=$temp[3];
		}
		else{
			@temp=split(/\s+/);
			$species=(split(/\./,$temp[1]))[0];
			push(@{$depth{$chrom}{$start}},$species);	
			@seq=split(//,$temp[6]);
			for $i(0..$#seq){
				$position=$start+$i;
				push(@{$seq{$chrom}{$position}{$species}},$seq[$i]);
				my %hash=();
				@{$seq{$chrom}{$position}{$species}}=sort grep{++$hash{$_}<2} @{$seq{$chrom}{$position}{$species}};
			}
		}
	}
	close(DATA);
	
	open(INPUT,">$work_dir/$target_chrom/stats/Hs.$chrom.$block.stats");
	for $chrom(sort keys %len){
		for $start(sort{$a<=>$b} keys %{$len{$chrom}}){
			$A0_len=$len{$chrom}{$start};
			my %hash=();
			my @array=@{$depth{$chrom}{$start}};
			@array=sort grep{++$hash{$_}<2} @array;
			
			for $i(0..$A0_len-1){
				$position=$start+$i;
				my @array_filter;
				$depth=$#array+1; ########stats speceis numbers in a block
				for $species(@array){
					if($#{$seq{$chrom}{$position}{$species}}==0 && ${$seq{$chrom}{$position}{$species}}[0] eq "-"){
						$depth=$depth-1;
						next;
					}
					push(@array_filter,$species);
					@array_filter=sort @array_filter;
				}
				
				$line=join(",",@array_filter);
				
				if(exists($phylop{$chrom}{$position})){
					print INPUT "$chrom\t$position\t$phylop{$chrom}{$position}\t$depth\t$line\n";
				}
				else{
					print INPUT "$chrom\t$position\tNull\t$depth\t$line\n";
					print "$chrom\t$position\n";
				}
			}
		}
	}
	close(INPUT);
}
