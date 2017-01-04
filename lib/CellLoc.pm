package CellLoc;

## cellular location information from Ing for human LLs

my $dataFileHome="./data.";
$dataFileHome=$ENV{'KRMLNM'/data} if (defined $ENV{'KRMLNM'});

my %llToLoc=();
my %locToLl=();

my $init=0;

sub _init
{
    return if ($init);
    $init=1;
    foreach my $file("$dataFileHome/centrosome.cat")
    {
	open(LOC,$file);
	while ($_ =<LOC>)
	{
	    chomp;
	    next if (/LOCUSLINK/);
	    my @lls=split(/\t/);
	    my $loc=shift(@lls); 
	    next if ($loc=~/human chromosome [0-9XY]/);
	    &addLocToLls($loc,@lls);
	}
    }
    open(IN,"/home/robison/categories/krmisc.cat");
    while ($_ = <IN>)
    {
	chomp;
	my @lls=split(/\t/);
	my $id=shift(@lls);
	my $cat=shift(@lls);
	if ($cat=~/co-immunoprecipitates with CENP-A chromatin/)
	{
	    $cat="CENP-A chromatin";
	}
	elsif ($id eq 'KRMisc81') 
	{
	    $cat="mitotic spindle";
	}
	elsif ($cat=~/14654843/)
	{
	    $cat="centrosome";
	}
	elsif ($cat=~/^(.*) localization in LIFEdb/)
	{
	    $cat=$1;
	}
	elsif ($cat=~/midbody/)
	{
	    $cat="midbody";
	}
	else
	{
	    next;
	}
	&addLocToLls($cat,@lls);
    }
    &_krTactical();
}
sub _krTactical
{
    foreach my $row(["centrosome",51199,22981,26524,10048,
		     1111,4957,4999,5347,5371,5901,7272,8195,8481,9859,10142,10564,10579,10636,10664,23265,23392,25932,55201,79848,164786,261734,345499,
		     219539,388403,51646,83719,29799,
		     4968,6103,57096,11200,55857,121441,10542,
		     23332,23122,4957,55697
],
		    ["midbody",11200],
		    ["nucleus",84962,4173,4175,9874],
		    ["cytoplasm",64506],
		    ["microtubule",10726],
		    ["Golgi",64689,10726,9088,10564,10142,23265],
		    ["anaphase promoting complex",51529,8881,996,8697,51433,51434,29882,29945,64682,25847,10393,991,246184,51343],
		    ['kinetochores',10403,11130,11243,11335,147841,23468,25936,57082,57405,79003,79980,83540,55143,147841,83540,79003,11243,57082,79172,
		     23332,23122
		    ],
		    ["spindle",79649,81610,55143,115106,10927,121441,
		     23332,23122
		    ],
		    ["nucleolus",219539,388403,51646,83719,29799]
   )
    {
	addLocToLls(@{$row})
    }
}
sub addLocToLls
{
    my @lls=@_;
    &_init() unless ($init);
    my $loc=shift(@lls);
    $loc="cytoplasm" if ($loc eq 'cytosol');
    foreach my $ll(@lls)
    {
	$llToLoc{$ll}={} unless (ref $llToLoc{$ll});
	$locToLl{$loc}={} unless (ref  $locToLl{$loc});
	$locToLl{$loc}->{$ll}=1;
	$llToLoc{$ll}->{$loc}=1;
    }
}
sub locations
{
    &_init() unless ($init);
    return keys %locToLl;
}
sub llForLocation
{
    my @locations=@_;
    &_init() unless ($init);
    my %ret=();
    foreach my $location(@locations)
    {
	if (ref $locToLl{$location})
	{
	    foreach my $ll(keys %{$locToLl{$location}})
	    {
		$ret{$ll}=1;
	    }
	}
    }
    return sort keys %ret;
}

sub locationsForLl
{
    my ($ll)=@_;
    &_init() unless ($init);
    my @ret=();
    @ret=keys %{$llToLoc{$ll}} if (ref $llToLoc{$ll});
    return @ret;
}

sub hasLocation
{
    my ($ll)=@_;
    &_init() unless ($init);
    return ref $llToLoc{$ll};
}
1;
