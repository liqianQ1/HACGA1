#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::Taxon;

###################################################################

my $opt_infile = shift;
my $opt_outprefix = shift;

###################################################################

my $opt_nodes = "/share/public/database/zhouy/Database/Taxonomy/nodes.dmp";
my $opt_names = "/share/public/database/zhouy/Database/Taxonomy/names.dmp";

my $unknown = "unknown";

###################################################################

sub main {

    my $query = &read_infile($opt_infile);
    my $dbh = &load_taxonomy($opt_nodes, $opt_names);

    my $qtax = &findout_tax($query, $dbh);
    my ($countDivision, $countTaxid, $mapCount) = &stat_taxonomy($qtax);

    my $outfile = "${opt_outprefix}.blast_nt.xls";
    &out_blast($qtax, $outfile);

    my $divisionFile = "${opt_outprefix}.stat_division.xls";
    my $taxidFile    = "${opt_outprefix}.stat_species.xls";
    &out_stat($countDivision, $countTaxid, $mapCount,
        $divisionFile, $taxidFile);

}

###################################################################

sub read_infile {
    my $infile = shift;

    my $query = {};

    ## dbj|AP014921.1|
    ## link , https://www.ncbi.nlm.nih.gov/nuccore/dbj%7CAP014921.1%7C
    open IN, "cut -f 1,2,13,16 $infile | " or die "error : can not open $infile\n";
    while (<IN>) {
        chomp;
        my ($qid, $sid, $taxid, $stitle) = split /\t/;

        die "error : repeat query id not allowed\n" if exists $query->{$qid};
        if ($taxid eq "N/A") {$taxid = $unknown};

        my $nt = $sid;
        $nt =~ s/\|/%7C/g;
        $nt = "https://www.ncbi.nlm.nih.gov/nuccore/$nt";

        $query->{$qid}->{sid} = $sid;
        $query->{$qid}->{nt} = $nt;
        $query->{$qid}->{tax} = $taxid;
        $query->{$qid}->{stitle} = $stitle;
    }
    close IN;

    return $query;
}

###################################################################

sub load_taxonomy {
    my ($nodes_db, $names_db) = @_;

    # Get one from a database
    my $dbh = Bio::DB::Taxonomy->new(-source   => 'flatfile',
                                 -directory=> '/tmp',
                                 -nodesfile=> $nodes_db,
                                 -namesfile=> $names_db);

    return $dbh;
}

###################################################################

sub findout_tax {
    my ($query, $dbh) = @_;

    my $qtax = {};

    foreach my $qid (keys %$query) {
        my $taxid = $query->{$qid}->{tax};

        my ($species, $division, $scientific_name);
        $species = $dbh->get_taxon(-taxonid => $taxid);

        if (defined $species && $species->division) {
            $division = $species->division;
            $scientific_name = $species->scientific_name;
        } else {
            $division = $unknown;
            $scientific_name = $unknown;
        }

        $qtax->{$qid} = {
            'sid'  =>  $query->{$qid}->{sid},
            'nt'  =>  $query->{$qid}->{nt},
            'taxid'  =>  $taxid,
            'division'  =>  $division,
            'taxname'  =>  $scientific_name,
            'stitle'  =>  $query->{$qid}->{stitle}
        };
    }

    return $qtax;
}

###################################################################

sub out_blast {
    my ($qtax, $outfile) = @_;

    open OUT,">",$outfile or die "error : can not write $outfile\n";
    foreach my $qid ( sort {
            $qtax->{$b}->{division} cmp $qtax->{$a}->{division} ||
            $qtax->{$b}->{taxid} cmp $qtax->{$a}->{taxid}
        } keys %$qtax ) {

        print OUT $qid;
        print OUT "\t".$qtax->{$qid}->{sid};
        print OUT "\t".$qtax->{$qid}->{nt};
        print OUT "\t".$qtax->{$qid}->{division};
        ## print OUT "\t".$qtax->{$qid}->{taxid};
        print OUT "\t".$qtax->{$qid}->{taxname};
        print OUT "\t".$qtax->{$qid}->{stitle};
        print OUT "\n";
    }
    close OUT;
}

###################################################################

sub stat_taxonomy {
    my ($qtax, $outstat) = @_;

    my $mapCount = scalar keys %$qtax;
    print STDERR "Total Mapped Reads : $mapCount\n";

    my $countDivision = {};
    my $countTaxid = {};
    foreach my $qid (keys %$qtax) { 
        my $division = $qtax->{$qid}->{division};
        $countDivision->{$division}++;

        my $taxid = $qtax->{$qid}->{taxname};
        $countTaxid->{$taxid}->{num}++;
        $countTaxid->{$taxid}->{division} = $division;
    }

    return ($countDivision, $countTaxid, $mapCount);
}

###################################################################

sub out_stat {
    my ($countDivision, $countTaxid, $mapCount, $divisionFile, $taxidFile) = @_;

    ###
    open OUT1,">",$divisionFile or die "error : can not write $divisionFile\n";
    print OUT1 "Division\tMapped_Count\tMapped_Percent\n";
    foreach my $division ( sort { $countDivision->{$b} <=> $countDivision->{$a}
        } keys %$countDivision) {

        print OUT1 $division;
        print OUT1 "\t".$countDivision->{$division};
        print OUT1 "\t".sprintf("%.2f", 100 * $countDivision->{$division} / $mapCount);
        print OUT1 "\n";
    }
    close OUT1;

    ###
    open OUT2,">",$taxidFile or die "error : can not write $taxidFile\n";
    print OUT2 "Species\tDivision\tMapped_Count\tMapped_Percent\n";
    foreach my $taxid ( 
            sort { $countTaxid->{$b}->{num} <=> $countTaxid->{$a}->{num}
        } keys %$countTaxid) {

        print OUT2 $taxid;
        print OUT2 "\t".$countTaxid->{$taxid}->{division};
        print OUT2 "\t".$countTaxid->{$taxid}->{num};
        print OUT2 "\t".sprintf("%.2f", 100 * $countTaxid->{$taxid}->{num} / $mapCount);
        print OUT2 "\n";
    }
    close OUT2;
}

###################################################################

main();


__END__


