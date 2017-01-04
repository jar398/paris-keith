use strict;
use IO::File;

package KrHomologene;

# most code really from Steve


my %geneToGroup=();
my %groups=();

my %TaxId2Species;
my %CommonName2Species;

my %Species = ();
foreach my $row([9606,"Homo sapiens","human"],
		[10090,'Mus musculus',"mouse"],
		[10116,'Rattus norvegicus','rat'],
		[8355,"Xenopus laevis","xenopus"],
		[4932,"Saccharomyces cerevisiae","yeast"],
		[6239,"Caenorhabditis elegans","worm"],
		[284812,"Schizosaccharomyces pombe","spombe"],
		[7227,"Drosophila melanogaster","fly"]
		# fly, pombe
	       )
{
    my ($taxId,$latin,$common)=@{$row};
    $Species{$latin}={'latin'=>$latin, 'common_name'=>$common,'tax_id'=>$taxId};
    $TaxId2Species{$taxId}=$Species{$latin};
    $CommonName2Species{$common}=$taxId;
}


my $initRun=0;


sub getHomologs
{
    my ($ll,$species)=@_;
    my $taxId=undef;
    if (defined $species)
    {
	if (defined $CommonName2Species{$species})
	{
	    $taxId=$CommonName2Species{$species} 
	}
	elsif (ref $Species{$species})
	{
	    $taxId=$Species{$species};
	}
    }
    &initHomologene() unless ($initRun);
    my @ret=();
    if ($ll>0 && ref $geneToGroup{$ll})
    {
	if (defined $taxId)
	{
	    if (ref $geneToGroup{$ll}->{'taxa'}->{$taxId})
	    {
		@ret=keys %{$geneToGroup{$ll}->{'taxa'}->{$taxId}};
	    }
	}
	else
	{
	    @ret=keys %{$geneToGroup{$ll}->{'genes'}};
	}
    }
    return @ret
}

sub initHomologene
{
    $initRun=1;
    &_parseHomologeneFile("$ENV{'DB'}/NCBI/HomoloGene/build52/homologene.data");
#"/genecatalog1/robison/txp/tgx/homologene.data.hum-rat"); # 
#    print STDERR scalar(keys %geneToGroup)," Homologene genes;\t",scalar(keys %hg)," groups","\n";
}

sub _parseHomologeneFile
{
    #    homologene.data is a tab delimited file containing the following
    #                    columns:
    #
    #             1) HID (HomoloGene group id)
    #             2) Taxonomy ID
    #             3) Gene ID
    #             4) Gene Symbol
    #             5) Protein gi
    #             6) Protein accession

    my($file)=@_;
#    print "Parsing Homologene file <$file>\n";
    my $fh = &_get_handle($file);



    while(<$fh>)
    {
        next if (/^\#/);
        chomp;
        my @f=split(/\t/,$_,-1);
        @f==6 or die();
        @f=map {s/^\s+//;s/\s+$//;$_} @f;
        #print "  line: ",join(",",@f),"\n";
        my($groupId,$tid,$geneId,$symbol,$prot_gi,$prot_acc)=@f;
        next unless (exists($TaxId2Species{$tid}));
        unless (defined($geneId) and $geneId=~/^\d+$/) {die();}
        $groups{$groupId}={'id'=>$groupId, 'genes'=>{},'taxa'=>{}} unless (ref $groups{$groupId});
	$groups{$groupId}->{'genes'}->{$geneId}=$tid;
	$groups{$groupId}->{'taxa'}->{$tid}={} unless (ref $groups{$groupId}->{'taxa'}->{$tid});
	$groups{$groupId}->{'taxa'}->{$tid}->{$geneId}=1;
	$geneToGroup{$geneId}=$groups{$groupId};
    }
    close($fh) or die();
    #print Dumper(\%Homologene);
    return 1;
}

sub _get_handle
{
    my($file)=@_;
    my $fh;
    -e $file or die("ERROR: file <$file> not found !\n");
    if ($file=~/\.gz$/){$fh=IO::File->new("gunzip -c $file|");}
    elsif ($file=~/\.Z$/){$fh=IO::File->new("uncompress -c $file|");}
    else {$fh=IO::File->new($file);}
    $fh or die("FATAL: unable to open file <$file> !\n");
    return $fh;
}

1;






