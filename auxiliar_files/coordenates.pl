###########################################
####### DEVELOPED by Lucas Miguel #########
### LGE - UNICAMP, Campinas, SP, Brazil ###
# contact: lucasmiguel@lge.ibi.unicamp.br #
###########################################

#!/usr/local/bin/perl

open FILE,"< $ARGV[0]";

my %cofactors=('phosphate'=>1,'H+'=>1,'CO2'=>1,'ATP'=>1,'ADP'=>1,'H'=>1,'NAD'=>1,'NADP'=>1,'NADH'=>1,'NADPH'=>1,'H2O'=>1);
while(<FILE>){
	my $line=$_;
	chomp $line;

	if($line=~/<AREA COORDS=\"(\d+),(\d+),(\d+),(\d+)\"/){
		my $id=$4;
		my $next=<FILE>;
		
		if($next=~/<b>Compound<\/b>: (.*)<br><b>Synonyms/){
			my $cmp=$1;
			if(not defined($cofactors{$cmp})){
				print $1."\t".$id."\n";
			}
		}
	}
}
			
close FILE;

