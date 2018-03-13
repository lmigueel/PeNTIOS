###########################################
####### DEVELOPED by Lucas Miguel #########
### LGE - UNICAMP, Campinas, SP, Brazil ###
# contact: lucasmiguel@lge.ibi.unicamp.br #
###########################################

#!/usr/bin/perl -w

use FindBin;
use File::Basename;

sub process_cmd {
        my ($cmd) = @_;

        my $sign = system($cmd);

        if ($sign) {
                die "Error: $cmd died with sign $ret";
        }

        return;
}
open LOG,"> output.log";


########################## DOWNLOAD AND PARSER PATHWAY FILE ################

my $dir = "$FindBin::RealBin";
my $time = localtime();
print LOG "[$time] Downloading $ARGV[0]....\n";

my $output = `grep \"$ARGV[0]\" $dir/sgd_pathway.file`;

($a,$b,$c)= split /\t/, $output;
chomp $a;
chomp $b;
chomp $c;

#auxiliar delay
system("wget '$c' -O test -o etest");
system("rm test");
system("rm etest");
###########


$cmd = "wget '$c' -O SAIDA -o exit";
&process_cmd($cmd);
$cmd = "rm exit";
&process_cmd($cmd);
$cmd="perl $dir/coordenates.pl 'SAIDA'  | sort -n -k 2 > $dir/coords.txt";
&process_cmd($cmd);

$time = localtime();
print LOG "[$time] Parsing expression values ... \n";
################# STORE NORMALIZED EXPRESSION VALUES ########

my $expression_file=$ARGV[1];

open EXP,"< $expression_file";

my %express;

while(<EXP>){
	my $line=$_;
        chomp $line;

	my @v=split(/\t/,$line);
        chomp $v[0];
        #print $v[0]."\t".$v[1]."\n";
        $express{$v[0]}=$v[1];
        
}
close EXP;

################################################################

$time = localtime();
print LOG "[$time] Parsing pathway file: $ARGV[0] ....\n";
	
open FILE,"< SAIDA";

my $pass=0;
my %genes;

while(<FILE>){
	my $line=$_;
	chomp $line;

	if($line=~/<MAP NAME=/){
		$pass=2;
	}
	if($pass==2){
		
		if($line=~/HREF=\"(.*)\"/){
			my $link=$1;
			if($link=~/RXN/){
                        	$cmd = "wget '$link' -O OUTPUT -o exit";
                        	&process_cmd($cmd);
                        	$cmd = "rm exit";
                        	&process_cmd($cmd);

				open ARQ,"< OUTPUT";
				my $flag=0;
				my $key="A";

				while(<ARQ>){	
					my $new=$_;
					chomp $new;
					if($new=~/Enzymes and Genes: /){
						$flag=1;
					}
					if($flag==1){
						if($new=~/<b>Gene:<\/b>(.*)<br> <b>Chromosome<\/b>/){
							my @v=split(/\s+/,$1);
							if($key eq "A"){
								$key=$v[2]."_".$v[1];
							}else{
								$key=$key.";".$v[2]."_".$v[1];
							}
							#print $key."\n";
						}
					}
					if($new=~/<b>Pathway<\/b>:/){#end genes
						$flag=0;
					}
				}
				close ARQ;										 
				$cmd="rm OUTPUT";
				&process_cmd($cmd);

				my $r=<FILE>;
				if($r=~/<b>Reaction<\/b>:(.*?)<br>/){
					my $react=$1;
					$genes{$key}=$react;
				}
                                elsif($r=~/<b>Reaction<\/b>:(.*?)\'/){
                                        my $react=$1;
                                        $genes{$key}=$react;
                                }

			}

		}
	}

				
}
close FILE;

$cmd="rm SAIDA";
&process_cmd($cmd);	

############ Aditional reactions ##################

my $new_r=$ARGV[2];

if($new_r eq "NO"){
}else{
	open REACT,"< $ARGV[2]";

	while (<REACT>){
		my $line=$_;
		chomp $line;

		my @list=split(/\t/,$line);
		my $id_gene="Y_".$list[0];
		my $reac=$list[1];
		my $expr=$list[2];

		$genes{$id_gene}=$reac;
		$express{$list[0]}=$expr;
		#print "*$id_gene*\t*$reac*\t*$expr*\n";
	}
	close REACT;
}
##################################### End PARSER ########################################

######################################### CREATE REACTIONS ##############################
$time = localtime();
print LOG "[$time] Creating reaction definitions ....\n";

## coordenates
open FILE,"< $dir/coords.txt";
my %coords;

while(<FILE>){
	my $line=$_;
	chomp $line;
	
	my @v=split(/\t/,$line);

	$coords{$v[0]}=$v[1];
}
close FILE;
###

#$cmd="rm $dir/coords.txt";
#&process_cmd($cmd);

my %hash_reaction; #reaction hash

foreach my $key (keys %genes){
	my $react=$genes{$key};
	#$react=~s/=/->/g;
	
	if($react=~/=/){
		my @r=split(/=/,$react);
                my @sub=split(/\+/,$r[0]);#substrates
                my @pro=split(/\+/,$r[1]);#products
		
		chomp $sub[0];
		$sub[0]=~s/\s+//g;
		chomp $pro[0];
		$pro[0]=~s/\s+//g;

		my $co1=$coords{$sub[0]};
		my $co2=$coords{$pro[0]};

		if($co1 > $co2){
			$react=~s/=/<-/g;
		}else{
			$react=~s/=/->/g;
		}
	}
	
	if($react=~/->/){

		my @r=split(/->/,$react);
                my @sub=split(/\+/,$r[0]);#substrates
                my @pro=split(/\+/,$r[1]);#products

                foreach my $a(@sub){	
			if($a eq " "){ next};
                        $a=~s/^\s+|\s+$//g;
			#print "Sub: $a\n";
			
			if(not defined($hash_reaction{$key})){
                        	$hash_reaction{$key}="s_".$a;
                        }else{
                                $hash_reaction{$key}=$hash_reaction{$key}.";s_".$a;
                        }
                }


                foreach my $a(@pro){
			if($a eq " "){ next};
			$a=~s/^\s+|\s+$//g;
			#print "Pro: $a\n";

                	$hash_reaction{$key}=$hash_reaction{$key}.";p_".$a;
                }
                                        
	
	}else{
                my @r=split(/<-/,$react);
                my @sub=split(/\+/,$r[1]);#substrates
                my @pro=split(/\+/,$r[0]);#products
                                
                foreach my $a(@sub){
			if($a eq " "){ next};

			$a=~s/^\s+|\s+$//g;
                        if(not defined($hash_reaction{$key})){
                                $hash_reaction{$key}="s_".$a;
                        }else{
                                $hash_reaction{$key}=$hash_reaction{$key}.";s_".$a;
                        }
                }       
                

                foreach my $a(@pro){
			if($a eq " "){ next};
			$a=~s/^\s+|\s+$//g;
                        $hash_reaction{$key}=$hash_reaction{$key}.";p_".$a;
                }

	}
}

##############################################################
$time = localtime();
print LOG "[$time] Normalizing data .... \n";

#################### NORMALIZED EXRESSION ##################

my $k=0; #count to array
my @array;
my %selected;#selected expression genes

foreach my $key (keys %hash_reaction){
	my @aux=split(/;/,$key);
	
	my $len=scalar @aux;
	my $gens;

	#log file with selected genes in case of multiple possible choice
	foreach my $a (@aux){
		my @v=split(/_/,$a);
		$a=$v[1];

		if($gens eq ""){
			$gens=$a;
		}else{
			$gens=$gens."_".$a;
		}

	}
	
	print LOG "The genes ($gens) have $len possibilities of choice ...\n";
	
	if($len ==1){

		#my @vec=split(/_/,$aux[0]);
		print LOG "Choosed was $gens with expression value equals to $express{$gens}\n";
		$selected{$key}=$express{$gens};
		$array[$k]=$express{$gens};
                $k++;

	}else{
		my $max=-999999;
		my $gen_max;

		foreach my $a (@aux){
			
			if(defined($express{$a})){
				if($express{$a} > $max){
					$max=$express{$a};
					$gen_max=$a;
				}
			}
		}

	        $array[$k]=$max;
	        $k++;
		print LOG "Choosed: gene $gen_max with highest expression value equalts to $max ... \n";	
		$selected{$key}=$max;
	}
}


@sort = sort{$a <=> $b}(@array);
my $len = scalar @sort;

if( $len % 2 == 0){
	my $x=($len/2)-1;
	my $num1=$sort[$x];
	my $num2=$sort[($len/2)];

        $sum=$num1+$num2;
        $med = $sum/2;
}else{
         $med=$sort[@sort/2];
}

foreach my $key (keys %selected){
        $selected{$key}=$selected{$key}/$med;
}

####################################################
$time = localtime();
print LOG "[$time] Generating SBML file ... \n";

print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<sbml xmlns=\"http://www.sbml.org/sbml/level2/version4\" level=\"2\" version=\"4\">
  <model>
    <listOfUnitDefinitions>
      <unitDefinition id=\"substance\">
        <listOfUnits>
          <unit kind=\"item\"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>
    <listOfCompartments>
      <compartment id=\"compartment\" spatialDimensions=\"0\"/>
    </listOfCompartments>
    <listOfSpecies>
";

my %list1; # hash of lists of species
foreach my $key (keys %hash_reaction){
        my @a=split(/;/,$hash_reaction{$key});

        foreach my $b (@a){
		#chomp $b;
                $b=~s/\/\///g;
                $b=~s/s_//g;
                $b=~s/p_//g;
                $b=~s/ /_/g;
                $b=~s/-/_/g;
                $b=~s/\)//g;
                $b=~s/\(//g;
                $b=~s/\[//g;
                $b=~s/,//g;

                if($b=~/^(\d+)_(.*)/){
                        $b="N".$1.$2;
                }

                if(not defined($list1{$b})){
                        print "     <species id=\"$b\" compartment=\"compartment\" initialAmount=\"0\" hasOnlySubstanceUnits=\"true\"/>\n";
                        $list1{$b}=1;
                }
        }
}

print "    </listOfSpecies>
    <listOfReactions>\n";

foreach my $key (keys %hash_reaction){
        my $ci;
        #print "gene=$key\t reaction=$hash_reaction{$key}\n";

        my @vec=split(/_/,$key);
        my $exp=$selected{$key};
	my @aux=split(/;/,$key);
	my $r_id="";

	foreach my $r (@aux){
		my @aux2=split(/_/,$r);
		if($r_id eq ""){
			$r_id=$aux2[1];
		}else{
			$r_id=$r_id."_".$aux2[1];
		}
	}
	#my $r_id=$key;
	#$r_id=~s/;/_/g;
	
	print "      <reaction id=\"$r_id\" reversible=\"false\">\n";
        print "       <listOfReactants>\n";
        my @a=split(/;/,$hash_reaction{$key});

        foreach my $b (@a){
                $b=~s/\/\///g;
                if($b=~/s_/){
                        $b=~s/s_//g;
                        $b=~s/ /_/g;
                        $b=~s/-/_/g;
                        $b=~s/\)//g;
                        $b=~s/\(//g;
                        $b=~s/\[//g;
                        $b=~s/,//g;

                        if($b=~/^(\d+)_(.*)/){
                                $b="N".$1.$2;
                        }
                        print "          <speciesReference species=\"$b\" stoichiometry=\"1\"/>\n";
                        if($ci eq ""){
                                $ci=$b;
                        }else{
                                $ci=$ci.";".$b; 
                        }
                }
        }
        print "       </listOfReactants>\n";
        print "        <listOfProducts>\n";

        foreach my $b (@a){
                $b=~s/\/\///g;
                if($b=~/p_/){
                        $b=~s/p_//g;
                        $b=~s/ /_/g;
                        $b=~s/-/_/g;
                        $b=~s/\)//g;
                        $b=~s/\(//g;
                        $b=~s/\[//g;
                        $b=~s/,//g;

                        if($b=~/^(\d+)_(.*)/){
                                $b="N".$1.$2;
                        }


                        print "          <speciesReference species=\"$b\" stoichiometry=\"1\"/>\n";
                }
        }
print"        </listOfProducts>
        <kineticLaw>
          <math xmlns=\"http://www.w3.org/1998/Math/MathML\">
            <apply>
              <times/>
";
print "              <cn>".$exp." </cn>\n";

        my @v=split(/;/,$ci);
        foreach my $x (@v){
                print "           <ci> $x </ci>\n";
        }
        print "            </apply>
          </math>
        </kineticLaw>
      </reaction>\n";

}
print "    </listOfReactions>
  </model>
</sbml>
";
$time = localtime();
print LOG "[$time] END with SUCCESS!\n";

close LOG;
	
