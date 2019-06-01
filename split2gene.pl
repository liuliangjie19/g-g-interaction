#!/usr/bin/perl -w
use strict;

sub getindex{
  (my $list1, my $list2) = @_;
  my @index = ();
  my $i=0;
  my $j=0;
  for($i=0;$i<=$#{$list1};$i++){
    for($j=0;$j<=$#{$list2};$j++){
      if (${$list1}[$i] eq ${$list2}[$j]){
        push(@index, $j);
        last;
      }
    }
  }
  return @index;
}

my $snpfile = $ARGV[0];
my $genofile = $ARGV[1];
#print $snpfile,$genofile;
my %gene2snp;
my $line;
my @linelist;
my @temp;

open(SNP, "<", $snpfile) || die "1 can not open the file!\n";

while(<SNP>){
  chomp;
  $line = $_;
  @linelist = split('\t', $line);
  if ($linelist[0] eq 'ID') {
    next;
  }
  if ($gene2snp{$linelist[2]}){
    push(@{$gene2snp{$linelist[2]}}, $linelist[5]);
  }
  else {
    @{$gene2snp{$linelist[2]}} = ($linelist[5]);
  }
}

close(SNP);

my $gene;
my @snplist;
my @indexlist;
open(GENO, "<", $genofile) || die "2 can not open the file!\n";
@linelist = <GENO>;
close(GENO);
$line = $linelist[0];
@snplist = split('\t', $line);
print @snplist;
foreach(keys(%gene2snp)){
  $gene = $_;
  @temp = @{$gene2snp{$gene}};
  $gene =~ s/\//_/g;
  open(OUT, ">", $gene.".txt") || die "3 can not open the file!\n";
  @indexlist = getindex(\@temp, \@snplist);
  print @indexlist,@temp;
  open(GENO, "<", $genofile) || die "4 can not open the file!\n";
  while(<GENO>){
    chomp;
    $line = $_;
    @linelist = split('\t', $line);
    #open(OUT, ">", $gene.".txt") || die "can not open the file!\n";
    print OUT "$linelist[0]\t$linelist[1]\t";
    foreach(@indexlist){
      print OUT "$linelist[$_]\t";
    }
    print OUT "\n";
  }
  close(OUT);
  close(GENO);
}
