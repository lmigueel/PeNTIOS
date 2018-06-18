###########################################
####### DEVELOPED by Lucas Miguel #########
### LGE - UNICAMP, Campinas, SP, Brazil ###
# contact: lucasmiguel@lge.ibi.unicamp.br #
###########################################

#!/usr/bin/env perl
use strict;
use warnings;

use FindBin;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

########################## PeNTIOS ##########################

 Requirements:

  --PathwayName <string>  : name of pathway in SGD between quotes (Ex: --Pathwayname "glycolysis")
	
  --expression_file <string> : file constaning egene and expression separated by <TAB>

  --output_file <string> : Output SBML file 

  --auxiliar_reactions <string> : New reactions to model (default: NO) 

  Example of usage:

  perl PeNTIOS.pl --PathwayName "glycolysis" --expression_file expression.txt --output_file output --auxiliar_reactions NO

  List all pathways in SGD:

  perl PeNTIOS.pl --pathway ALL

  List all pathways in SGD containing some word (ex: glucose):
 
  perl PeNTIOS.pl --pathway glucose
  
###############################################################

_EOUSAGE_

        ;

my ($sgd_id,$expression_file,$output_file);
my $path;
my $auxiliar;

&GetOptions(

    ## general opts
    "PathwayName=s" => \$sgd_id,
    "output_file=s" => \$output_file,
    "expression_file=s" => \$expression_file,
    "pathway=s" => \$path,
    "auxiliar_reactions=s" => \$auxiliar,

);

my $dir = "$FindBin::RealBin";

# Verify options:

if($path){

	if($path eq "ALL"){
		system("cut -f 2 $dir/auxiliar_files/sgd_pathway.file");
	}else{
		system("grep \"$path\" $dir/auxiliar_files/sgd_pathway.file | cut -f 2");
	}
}

unless ($sgd_id && $output_file && $expression_file) {
        die $usage;
}

main: {
	my $cmd;

	if($auxiliar ne "NO"){
		$cmd="perl $dir/auxiliar_files/generate_file.pl \"$sgd_id\" $expression_file $auxiliar > $output_file.sbml";
		&process_cmd($cmd);
	}else{
		$cmd="perl $dir/auxiliar_files/generate_file.pl \"$sgd_id\" $expression_file NO > $output_file.sbml";
		&process_cmd($cmd);
	}

exit(0);

}


sub process_cmd {
        my ($cmd) = @_;
        my $sign = system($cmd);

        if ($sign) {
                die "Error: $cmd died with sign $sign";
        }

        return;
}

